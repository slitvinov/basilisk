/**
# A vertically-Lagrangian multilayer solver for free-surface flows

This is the implementation of the solver described in [Popinet,
2020](/Bibliography#popinet2020). The system of $n$ layers of
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

#define BGHOSTS 2
#define LAYERS 1
#include "run.h"
#include "bcg.h"

scalar zb[], eta, h;
vector u;
double G = 1., dry = 1e-12 /* fixme: 1e-4 before */, CFL_H = nodata;
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

#include "run.h"

/**
Other fields, such as $\mathbf{u}_k$ here, are added by this
event. */

event defaults (i = 0)
{

  /**
  The (velocity) CFL is limited by the unsplit advection scheme, so is
  dependent on the dimension. The (gravity wave) CFL is set to 1/2 (if
  not already set by the user). */
  
  CFL = 1./(2.*dimension);
  if (CFL_H == nodata)
    CFL_H = 0.5;  
  
  u = new vector[nl];

  // fixme: hu, ha and hf are used only by diffusion.h
  
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

  /**
  We setup the default display. */

  display ("squares (color = 'h > 0 ? eta : nodata', spread = -1);");
}

/**
After user initialisation, we define the free-surface height $\eta$
and initial (face) acceleration. */

event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
  }
  boundary (all);
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

static void advect_list (double dt, scalar * tracers)
{
#if !NOCFL // fixme: does not work for implicit-ml.tst
  foreach_face()
    foreach_layer() {
      if (hu.x[]*dt/(Delta*cm[-1]) > CFL*h[-1])
	hu.x[] = CFL*h[-1]*Delta*cm[-1]/dt;
      else if (- hu.x[]*dt/(Delta*cm[]) > CFL*h[])
	hu.x[] = - CFL*h[]*Delta*cm[]/dt;
    }
  boundary ((scalar *){hu});
#endif
  
  /**
  ## Advection */

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

#if 0 // dimension > 1
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
    h_i^{n+1} & = h_i^n + \Delta t \frac{(hu)_{i+1/2} - (hu)_{i-1/2}}{\Delta},\\
    s_i^{n+1} & = \frac{(hs)^\star_i}{h_i^{n+1}}
    \end{aligned}
    $$
    */

    foreach() {
      double h1 = h[];
      foreach_dimension()
	h1 += dt*(hu.x[] - hu.x[1])/(Delta*cm[]);
      assert (h1 >= 0.);
      if (h1 < dry) {
	for (scalar f in tracers)
	  f[] = 0.;
      }
      else
	for (scalar f in tracers)
	  f[] /= h1;
    }
  }
  boundary (tracers);
}

void advect_h (double dt)
{
  foreach()
    foreach_layer() {
      double h1 = h[];
      foreach_dimension()
	h1 += dt*(hu.x[] - hu.x[1])/(Delta*cm[]);
      //      assert (h1 >= 0.);      
      h[] = max(h1,0.);
    }
  boundary ({h});
}

trace
void advect (double dt)
{
  advect_list (dt, tracers);
  advect_h (dt);
}

#define STABILITY 1

#if STABILITY
double dtmax;

event set_dtmax (i++,last) dtmax = DT;

static bool hydrostatic = true;

event stability (i++,last)
{
  foreach (reduction (min:dtmax)) {
    double H = 0.;
    foreach_layer()
      H += h[];
    if (H > dry) {
      double cp = hydrostatic ? sqrt(G*H) : sqrt(G*Delta*tanh(H/Delta));
      cp /= CFL_H;
      foreach_layer()
	foreach_dimension() {
	  double c = cp + fabs(u.x[]/CFL);
	  if (c > 0.) {
	    double dt = 2.*cm[]/(fm.x[] + fm.x[1])*Delta/c;
	    if (dt < dtmax)
	      dtmax = dt;
	  }
        }
    }
  }
  dt = dtnext (dtmax);
}
#endif // STABILITY

#define gmetric(i) (2.*fm.x[i]/(cm[i] + cm[i-1]))

event face_fields (i++, last)
{
  double dtmax = DT;
  foreach_face (reduction (min:dtmax)) {
    double ax = - gmetric(0)*G*(eta[] - eta[-1])/Delta;
    double H = 0., um = 0.;
    foreach_layer() {
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      hu.x[] = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;

#if 1 // important for stability of gaussian-hydro.tst
      // also has a significant influence on the amplitude of
      // the non-hydrostatic waves for gaussian.tst
      // setting a larger prefactor (e.g. 2.*dt) also stabilises
      double hff;
      if (h[] + h[-1] < dry) // fixme: useful??
	hff = 0.;
      else {
	// fixme: dt needs to be multiplied by (1. - theta_H)
	//	double un = (1. - theta_H)*dt*(hu.x[] + (1. - theta_H)*dt*ax)/Delta,
	//	  a = sign(un);
	double un = dt*(hu.x[] + dt*ax)/Delta, a = sign(un);
	int i = - (a + 1.)/2.;
	double g = h.gradient ? h.gradient (h[i-1], h[i], h[i+1])/Delta :
	  (h[i+1] - h[i-1])/(2.*Delta);
	hff = h[i] + a*(1. - a*un)*g*Delta/2.;
      }
#endif
      hf.x[] = fm.x[]*hff;

#if 0
      hu.x[] *= hf.x[];
      if (hu.x[]/cm[-1] > CFL*um*h[-1])
	um = hu.x[]/(CFL*h[-1]*cm[-1]);
      else if (- hu.x[]/cm[] > CFL*um*h[])
	um = - hu.x[]/(CFL*h[]*cm[]);
#else
      if (fabs(hu.x[]) > um)
	um = fabs(hu.x[]);
      hu.x[] *= hf.x[];
#endif
      ha.x[] = hf.x[]*ax;
 
      H += hff;
    }
    if (H > 0.) {
      double c = um/CFL + hydrostatic ? sqrt(G*H)/CFL_H :
	sqrt(G*Delta*tanh(H/Delta))/CFL_H;
      double dt = min(cm[], cm[-1])*Delta/(c*fm.x[]);
      if (dt < dtmax)
	dtmax = dt;
    }
  }
  boundary ((scalar *){ha, hu, hf});
#if !STABILITY  
#if 1
  static double previous = 0.;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
#endif
  dt = dtnext (dtmax);
#endif // !STABILITY
}

event coriolis (i++, last);
 
event pre_acceleration (i++, last);

event acceleration (i++, last)
{
  foreach_face()
    foreach_layer() 
      hu.x[] += dt*ha.x[];
  boundary ((scalar *){hu});
  
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
  boundary ((scalar *){u});

#if 1
  //  delete ((scalar *){ha});
  advect (dt);
  //  delete ((scalar *){hu, hf});
#endif
}

scalar deta[]; // fixme: just for diagnostic

event advection (i++, last)
{
  
  /**
  Finally the free-surface height $\eta$ is updated and the boundary
  conditions are applied. */

  foreach() {
    double etap = zb[];
    foreach_layer()
      etap += h[];
    deta[] = etap - eta[];
    eta[] = etap;
  }
  
  scalar * list = list_copy ({h, eta});
  for (scalar s in tracers)
    list = list_append (list, s);
  boundary (list);
  free (list);
}

/**
## Cleanup

The fields and lists allocated in [`defaults()`](#defaults0) above
must be freed at the end of the run. */
   
event cleanup (t = end, last)
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
#if 0  
  return H > dry ? sqrt(G/H)*(zb[] + H - ref) : 0.;
#else
  return sqrt (G*max(H,0.)) - sqrt(G*max(ref - zb[], 0.));
#endif
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

#if dimension == 2

/**
# Fluxes through sections

These functions are typically used to compute fluxes (i.e. flow rates)
through cross-sections defined by two endpoints (i.e. segments). Note
that the orientation of the segment is taken into account when
computing the flux i.e the positive normal direction to the segment is
to the "left" when looking from the start to the end. 

This can be expressed mathematically as:
$$
\text{flux}[k] = \int_A^B h_k\mathbf{u}_k\cdot\mathbf{n}dl
$$
with $A$ and $B$ the endpoints of the segment, $k$ the list index
(typically the layer index), $\mathbf{n}$ the oriented segment unit
normal and $dl$ the elementary length. The function returns the sum
(over $k$) of all the fluxes. */

double segment_flux (coord segment[2], double * flux, scalar h, vector u)
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  for (int l = 0; l < nl; l++)
    flux[l] = 0.;
  foreach_segment (segment, p) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);    
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      foreach_layer()
	flux[_layer] += dl/2.*
	interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.})*
	(m.x*interpolate_linear (point, (struct _interpolate)
				 {u.x, a.x, a.y, 0.}) +
	 m.y*interpolate_linear (point, (struct _interpolate)
				 {u.y, a.x, a.y, 0.}));
    }
  }
  // reduction
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, flux, nl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  double tot = 0.;
  for (int l = 0; l < nl; l++)
    tot += flux[l];
  return tot;
}

/**
A NULL-terminated array of *Flux* structures passed to
*output_fluxes()* will create a file (called *name*) for each
flux. Each time *output_fluxes()* is called a line will be appended to
the file. The line contains the time and the value of the flux for
each $h$, $u$ pair in the lists. The *desc* field can be filled with a
longer description of the flux. */

typedef struct {
  char * name;
  coord s[2];
  char * desc;
  FILE * fp;
} Flux;

void output_fluxes (Flux * fluxes, scalar h, vector u)
{
  for (Flux * f = fluxes; f->name; f++) {
    double flux[nl];
    double tot = segment_flux (f->s, flux, h, u);
    if (pid() == 0) {
      if (!f->fp) {
	f->fp = fopen (f->name, "w");
	if (f->desc)
	  fprintf (f->fp, "%s\n", f->desc);
      }
      fprintf (f->fp, "%g %g", t, tot);
      for (int i = 0; i < nl; i++)
	fprintf (f->fp, " %g", flux[i]);
      fputc ('\n', f->fp);
      fflush (f->fp);
    }
  }
}

#endif // dimension == 2
