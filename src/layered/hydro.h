/**
# A vertically-Lagrangian multilayer solver for free-surface flows

This is the implementation of the solver described in [Popinet,
2019](/Bibliography#popinet2019). The system of $n$ layers of
incompressible fluid pictured below is modelled as the set of
(hydrostatic) equations:
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
\end{aligned}
$$
with $\mathbf{u}_k$ the velocity vector, $h_k$ the layer thickness,
$g$ the acceleration of gravity and $\eta$ the free-surface
height. The non-hydrostatic pressure $\phi_k$ and vertical velocity
$w_k$ can be added using [the non-hydrostatic extension](nh.h).

![Definition of the $n$ layers.](../figures/layers.svg){ width="60%" }

## Fields and parameters

The `zb` and `eta` fields define the bathymetry and free-surface
height. The `h` and `u` fields define the layer thicknesses and velocities.

The acceleration of gravity is `G` and `dry` controls the minimum
layer thickness. The hydrostatic CFL criterion is defined by `CFL_H`.

The `gradient` pointer gives the function used for limiting.

The `hu`, `ha` and `hf` fields define the face flux $hu$,
height-weighted face acceleration and face height for each layer.

`tracers` is a list of tracers for each layer. By default it contains
only the components of velocity. */

#include "run.h"
#include "bcg.h"

scalar zb[], eta, h;
vector u;
double G = 1., dry = 1e-4, CFL_H = 1.;
double (* gradient) (double, double, double) = minmod2;
bool linearised = false;

vector hu, ha, hf;
scalar * tracers = NULL;

/**
## Setup

We first define refinement and restriction functions for the
free-surface height `eta` which is just the sum of all layer
thicknesses and of bathymetry. */

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child() {
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
  }
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[];
  foreach_layer()
    eta[] += h[];
}
#endif // TREE

/**
The allocation of fields for each layer is split in two parts, because
we need to make sure that the layer thicknesses and $\eta$ are
allocated first (in the case where other layer fields are added, for
example for the [non-hydrostatic extension](nh.h)). */

event defaults0 (i = 0)
{
  assert (nl > 0);
  h = new scalar[nl];
  h.gradient = gradient;
#if TREE
  h.refine = h.prolongation = refine_linear;
  h.restriction = restriction_volume_average;
#endif
  eta = new scalar;
  reset ({h, zb}, 0.);

  /**
  We set the proper gradient and refinement/restriction functions. */
  
  zb.gradient = gradient;
  eta.gradient = gradient;
#if TREE
  zb.refine = zb.prolongation = refine_linear;
  zb.restriction = restriction_volume_average;
  eta.prolongation = refine_linear;
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
#endif // TREE
}

/**
Other fields, such as $\mathbf{u}_k$ here, are added by this
event. */

event defaults (i = 0)
{
  u = new vector[nl];
  hu = new face vector[nl];
  ha = new face vector[nl];
  hf = new face vector[nl];
  reset ({u, hu, ha, hf}, 0.);
    
  if (!linearised)
    foreach_dimension()
      tracers = list_append (tracers, u.x);
  
  /**
  The gradient and prolongation/restriction functions are set for all
  tracer fields. */

  for (scalar s in tracers) {
    s.gradient = gradient;
#if TREE
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
#endif
  }
}

/**
After user initialisation, we define the free-surface height $\eta$
and initial acceleration. */

double dtmax;

event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
  }
  boundary (all);

  dt = 0;
  event ("acceleration");
  dtmax = DT;
  event ("stability");
  event ("acceleration");
}

/**
## Stability condition

The maximum timestep is set using the standard (hydrostatic) CFL condition
$$
\Delta t < \text{CFL}\frac{\Delta}{|\mathbf{u}|+\sqrt{gH}}
$$
or the modified (non-hydrostatic) condition ([Popinet,
2019](/Bibliography#popinet2019))
$$
\Delta t < 
\text{CFL}\frac{\Delta}{|\mathbf{u}|+\sqrt{g\Delta H\text{tanh}(H/\Delta)}}
$$
*/

event set_dtmax (i++,last) dtmax = DT;

static bool hydrostatic = true;

event stability (i++,last)
{
  foreach_face (reduction (min:dtmax)) {
    double Hf = 0.;
    foreach_layer()
      Hf += hf.x[];
    if (Hf > dry) {
      Hf /= fm.x[];
      double cp = hydrostatic ? sqrt(G*Hf) : sqrt(G*Delta*tanh(Hf/Delta));
      cp /= CFL_H;
      foreach_layer() {
	double c = hf.x[]*cp + fabs(hu.x[]);
	if (c > 0.) {
	  double dt = hf.x[]/fm.x[]*cm[]*Delta/c;
	  if (dt < dtmax)
	    dtmax = dt;
	}
      }
    }
  }
  dt = dtnext (CFL*dtmax);
}

/**
## Advection (and diffusion)

This section first solves
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = 0, \\
  \partial_t \left( h s \right)_k + \mathbf{{\nabla}} \cdot \left(
  h s  \mathbf{u} \right)_k & = 0,
\end{aligned}
$$
where $s_k$ is any other tracer field added to the `tracers` list. */

event advection_term (i++,last)
{
  if (grid->n == 1) // to optimise 1D-z models
    return 0;
  
  face vector flux[];
  
  foreach_layer() {
    
    /**
    We compute the flux $(shu)_{i+1/2,k}$ for each tracer $s$, using a
    variant of the BCG scheme. */
    
    for (scalar s in tracers) {
      foreach_face() {
	double un = dt*hu.x[]/(hf.x[]*Delta + dry), a = sign(un);
	int i = -(a + 1.)/2.;
	double g = s.gradient ?
	  s.gradient (s[i-1], s[i], s[i+1])/Delta :
	  (s[i+1] - s[i-1])/(2.*Delta);
	double s2 = s[i] + a*(1. - a*un)*g*Delta/2.;
#if dimension > 1
	if (hf.y[i] + hf.y[i,1] > dry) {
	  double vn = (hu.y[i] + hu.y[i,1])/(hf.y[i] + hf.y[i,1]);
	  double syy = (s.gradient ? s.gradient (s[i,-1], s[i], s[i,1]) :
			vn < 0. ? s[i,1] - s[i] : s[i] - s[i,-1]);
	  s2 -= dt*vn*syy/(2.*Delta);
	}
#endif
	flux.x[] = s2*hu.x[];
      }
      boundary_flux ({flux});
      
      /**
      We compute $(hs)^\star_i = (hs)^n_i + \Delta t 
      [(shu)_{i+1/2} -(shu)_{i-1/2}]/\Delta$. */
      
      foreach() {
	s[] *= h[];
	foreach_dimension()
	  s[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
      }
    }

    /**
    We then obtain $h^{n+1}$ and $s^{n+1}$ using
    $$
    \begin{aligned}
    h_i^{n+1} & = h_i^n + \Delta t \frac{(hu)_{i+1/2} - (hu)_{i-1/2}}{\Delta}, \\
    s_i^{n+1} & = \frac{(hs)^\star_i}{h_i^{n+1}}
    \end{aligned}
    $$
    */

    foreach() {
      foreach_dimension()
	h[] += dt*(hu.x[] - hu.x[1])/(Delta*cm[]);
      if (h[] < dry) {
	for (scalar f in tracers)
	  f[] = 0.;
      }
      else
	for (scalar f in tracers)
	  f[] /= h[];
    }
  }
  
  /**
  Finally the free-surface height $\eta$ is updated and the boundary
  conditions are applied. */

  foreach() {
    double H = 0.;
    foreach_layer()
      H += h[];
    eta[] = zb[] + H;
  }
  
  scalar * list = list_copy ({h, eta});
  for (scalar s in tracers)
    list = list_append (list, s);
  boundary (list);
  free (list);
}

/**
[Vertical remapping](remap.h) is applied here if necessary. */

event remap (i++, last);

/**
Vertical diffusion (including viscosity) is added by this code. */

#include "diffusion.h"

/**
## Mesh adaptation

Plugs itself here. */

#if TREE
event adapt (i++,last);
#endif

/**
## Acceleration */

event acceleration (i++,last)
{
  trash ({hu, ha, hf});
  foreach_face() {
    
    /**
    We compute the acceleration as
    $$
    a_{i + 1 / 2, k} = - g \frac{\eta_{i + 1, k} - \eta_{i, k}}{\Delta}
    $$
    taking into account [metric
    terms](/src/README#general-orthogonal-coordinates).*/
    
    double ax = - fm.x[]*G*(eta[] - eta[-1])/((cm[] + cm[-1])*Delta/2.);
    foreach_layer() {

      /**
      The face velocity is first computed as
      $$
      u_{i + 1 / 2, k} = \frac{(hu)_{i + 1, k} + 
      (hu)_{i, k}}{h_{i + 1, k} + h_{i, k}} + \Delta ta_{i + 1 / 2, k}
      $$
      with checks to avoid division by zero. */

      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      hu.x[] = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;

      /**
      The face height `hf` is reconstructed using Bell-Collela-Glaz-style
      upwinding. */
	
      double un = dt*(hu.x[] + dt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = h.gradient ? h.gradient (h[i-1], h[i], h[i+1])/Delta :
	(h[i+1] - h[i-1])/(2.*Delta);
      hf.x[] = h[i] + a*(1. - a*un)*g*Delta/2.;

      /**
      The height-weighted face acceleration $ha_{i+1/2}$ and
      face flux $hu_{i+1/2}$ are computed as */

      if (hf.x[] < dry)
	ha.x[] = hu.x[] = hf.x[] = 0.;
      else {
	hf.x[] *= fm.x[];
	ha.x[] = hf.x[]*ax;	
	hu.x[] = hf.x[]*(hu.x[] + dt*ax);
      }
    }
  }
  boundary ((scalar *){hu, ha, hf});
}

/**
Finally the acceleration is used to obtain the velocity field at time
$n+1$ as
$$
u^{n + 1}_{i, k} \leftarrow u^{n + 1}_{i, k} + \Delta t \frac{ha_{i + 1 / 2,
k} + ha_{i - 1 / 2, k}}{h_{i + 1 / 2,k} + h_{i - 1 / 2, k}}
$$
The missing shallow-water metric terms $(f_G v, -f_G u)$, with
$$
f_G \equiv \frac{v \partial_{\lambda} m_{\theta} - u \partial_{\theta}
m_{\lambda}}{m_{\lambda} m_{\theta}}
$$
are also added here ([Popinet, 2011](/src/references.bib#popinet2011)). */

event pressure (i++,last)
{
  foreach()
    foreach_layer() {
      foreach_dimension()
        u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
#if dimension == 2
      // metric terms
      double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
      double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
      double ux = u.x[], uy = u.y[];
      double fG = uy*dmdl - ux*dmdt;
      u.x[] += dt*fG*uy;
      u.y[] -= dt*fG*ux;
#endif // dimension == 2
    }
  boundary ((scalar *) {u});
}

/**
## Cleanup

The fields and lists allocated in [`defaults()`](#defaults0) above
must be freed at the end of the run. */
   
event cleanup (i = end, last)
{
  delete ({eta, h, u, hu, ha, hf});
  free (tracers), tracers = NULL;
}

/**
# Hydrostatic vertical velocity

For the hydrostatic solver, the vertical velocity is not defined by
default since it is usually not required. The function below can be
applied to compute it using the (diagnostic) incompressibility condition
$$
\mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k + \left[ w -
\mathbf{u} \cdot \mathbf{{\nabla}} (z) \right]_k = 0
$$
*/

void vertical_velocity (scalar w)
{
  foreach() {
    double dz = zb[1] - zb[-1];
    double wm = 0.;
    foreach_layer() {
      w[] = wm + (hu.x[] + hu.x[1])/(hf.x[] + hf.x[1] + dry)*
	(dz + h[1] - h[-1])/(2.*Delta);
      if (point.l > 0)
	foreach_dimension()
	  w[] -= (hu.x[0,0,-1] + hu.x[1,0,-1])
	  /(hf.x[0,0,-1] + hf.x[1,0,-1] + dry)*dz/(2.*Delta);
      foreach_dimension()
	w[] -= (hu.x[1] - hu.x[])/Delta;
      dz += h[1] - h[-1], wm = w[];
    }
  }
  boundary ({w});
}

/**
# "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (*ref*). */

double _radiation (Point point, double ref, scalar s)
{
  double H = 0.;
  foreach_layer()
    H += h[];
  return H > dry ? sqrt(G/H)*(zb[] + H - ref) : 0.;
}
  
#define radiation(ref) _radiation(point, ref, _s)

/**
# Conservation of water surface elevation

We re-use some generic functions. */

#include "elevation.h"

/**
But we need to re-define the water depth refinement function. */

#if TREE
static void refine_layered_elevation (Point point, scalar h)
{

  /**
  We first check whether we are dealing with "shallow cells". */
  
  bool shallow = zb[] > default_sea_level;
  foreach_child()
    if (zb[] > default_sea_level)
      shallow = true, break;

  /**
  If we do, refined cells are just set to the default sea level. */
  
  if (shallow)
    foreach_child()
      h[] = max(0., default_sea_level - zb[]);

  /**
  Otherwise, we use the surface elevation of the parent cells to
  reconstruct the water depth of the children cells. */
  
  else {
    double eta = zb[] + h[];   // water surface elevation  
    coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
	g.x = gradient (zb[-1] + h[-1], eta, zb[1] + h[1])/4.;
    else
      foreach_dimension()
	g.x = (zb[1] + h[1] - zb[-1] - h[-1])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta;
      foreach_dimension()
	etac += g.x*child.x;
      h[] = max(0., etac - zb[]);
    }  
  }
}

/**
We overload the `conserve_elevation()` function. */

void conserve_layered_elevation (void)
{
  h.refine  = refine_layered_elevation;
  h.prolongation = prolongation_elevation;
  h.restriction = restriction_elevation;
  boundary ({h});
}

#define conserve_elevation() conserve_layered_elevation()

#endif // TREE

#include "gauges.h"
