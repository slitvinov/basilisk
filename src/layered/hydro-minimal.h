/**
# Minimal hydro (for testing 1D models only)

This is only temporary and will be removed when the layered
formulation has been rewritten. */

#include "run.h"
#include "bcg.h"

int nl = 1;
attribute { int l; }
scalar zb[], eta;
scalar * hl = NULL, ** tracers = NULL;
vector * ul = NULL;
double G = 1., dry = 1e-6, CFL_H = 1.;
double (* gradient) (double, double, double) = minmod2;
bool linearised = false;

vector * ufl = NULL, * al = NULL;
scalar eta_a;

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
    for (scalar h in hl)
      eta[] += h[];
  }
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[];
  for (scalar h in hl)
    eta[] += h[];
}
#endif // TREE

/**
The allocation of fields for each layer is split in two parts, because
we need to make sure that the layer thicknesses and $\eta$ are
allocated first (in the case where other layer fields are added, for
example for the [non-hydrostatic extension](nh.h)). */

static scalar new_layered_scalar (const char * name, int l)
{
  scalar s = new scalar;
  char lname[80];
  sprintf (lname, "%s%d", name, l);
  free (s.name);
  s.name = strdup (lname);
  s.l = l;
  return s;
}

static vector new_layered_vector (const char * name, int l)
{
  vector v = new vector;
  foreach_dimension() {
    struct { char * x, * y; } c = {"x", "y"};
    char lname[80];
    sprintf (lname, "%s%d.%s", name, l, c.x);
    free (v.x.name);
    v.x.name = strdup (lname);
    v.x.l = l;
  }
  return v;
}

event defaults0 (i = 0)
{
  assert (hl == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar h = new_layered_scalar ("h", l);
    h.gradient = gradient;
#if TREE    
    h.refine = h.prolongation = refine_linear;
    h.restriction = restriction_volume_average;
#endif
    hl = list_append (hl, h);
  }
  eta = new scalar;
  eta_a = eta;
  reset (hl, 0.);
  reset ({zb}, 0.);

#if 0
  eta[left] = dirichlet ((4.*eta[] - 3.*eta[1] + eta[2])/2.);
  //  eta[left] = 3.*(eta[] - eta[1]) + eta[2];
  eta[right] = dirichlet ((4.*eta[] - 3.*eta[-1] + eta[-2])/2.);
  //  eta[right] = 3.*(eta[] - eta[-1]) + eta[-2];
#if dimension > 1 // fixme: foreach_dimension() does not work
  eta[bottom] = dirichlet ((4.*eta[] - 3.*eta[0,1] + eta[0,2])/2.);
  //  eta[bottom] = 3.*(eta[] - eta[0,1]) + eta[0,2];
  eta[top] = dirichlet ((4.*eta[] - 3.*eta[0,-1] + eta[0,-2])/2.);
  //  eta[top] = 3.*(eta[] - eta[0,-1]) + eta[0,-2];
#endif
#endif

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

  assert (!tracers);
  tracers = calloc (nl, sizeof(scalar *));
}

/**
Other fields, such as $\mathbf{u}_k$ here, are added by this
event. */

event defaults (i = 0)
{
  assert (ul == NULL && ufl == NULL && al == NULL);
  for (int l = 0; l < nl; l++) {
    vector u = new_layered_vector ("u", l);
    vector uf = new face vector;
#if 0
    uf.n[left] = 0.;
    uf.n[right] = 0.;
#if dimension > 1 // fixme: foreach_dimension() does not work
    uf.n[top] = 0.;
    uf.n[bottom] = 0.;
#endif
#endif

    vector a = new face vector;
    ul = vectors_append (ul, u);
    ufl = vectors_append (ufl, uf);
    al = vectors_append (al, a);
    foreach_dimension() {
      uf.x.l = a.x.l = l;
      if (!linearised)
	tracers[l] = list_append (tracers[l], u.x);
    }
  }
  reset (ul, 0.);
  reset (ufl, 0.);
  reset (al, 0.);

  /**
  The gradient and prolongation/restriction functions are set for all
  tracer fields. */

  for (int l = 0; l < nl; l++)
    for (scalar s in tracers[l]) {
      s.gradient = gradient;
#if TREE
      s.refine = s.prolongation = refine_linear;
      s.restriction = restriction_volume_average;
#endif
    }
}

/**
After user initialisation, we define the initial face velocities `uf`
and free-surface height $\eta$. */

double dtmax;

#if 1
#define face_is_wet() ((H > dry && Hm > dry) ||				\
		       (H > dry && eta[] >= zb[-1]) ||			\
		       (Hm > dry && eta[-1] >= zb[]))
#else
#define face_is_wet() ((H > dry && Hm > dry) ||				\
		       (H > dry && eta[] > (zb[] + zb[-1])/2.) ||	\
		       (Hm > dry && eta[-1] > (zb[] + zb[-1])/2.))
#endif

event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    for (scalar h in hl)
      eta[] += h[];
  }
  boundary (all);
  
  trash ((scalar *)ufl);
  foreach_face() {
    double H = 0., Hm = 0.;
    for (scalar h in hl)
      H += h[], Hm += h[-1];
    scalar h;
    face vector uf;
    vector u;
    if (face_is_wet()) { // fixme: unify this
      double H = 0., Hm = 0.;
      for (h,uf,u in hl,ufl,ul) {
	H += h[], Hm += h[-1];
	if ((H > dry && Hm > dry) ||
	    (H > dry && zb[] + H >= zb[-1]) ||
	    (Hm > dry && zb[-1] + Hm >= zb[]))
	  uf.x[] = fm.x[]*(max(h[],0.)*u.x[] +
			   max(h[-1],0.)*u.x[-1])/(max(h[],0.) +
						   max(h[-1],0.) + dry);
	else
	  uf.x[] = 0.;
      }
    }
    else
      for (vector uf in ufl)
	uf.x[] = 0.;
  }
  boundary ((scalar *)ufl);
 
  dtmax = DT;
  event ("stability");
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
    double H = 0., Hm = 0.;
#if !TREE
    for (scalar h in hl)
      H += h[], Hm += h[-1];
#else // TREE
    scalar h;
    vector uf;
    for (h,uf in hl,ufl) {
      H += h[], Hm += h[-1];
    
      /**
      We check and enforce that fluxes between dry and wet cells are
      zero. This is necessary only when using adaptive refinement,
      otherwise this is guaranteed by the way face velocities are
      computed in the [acceleration section](#acceleration). */
      
      if ((h[] < dry && uf.x[] < 0.) ||
	  (h[-1] < dry && uf.x[] > 0.))
	uf.x[] = 0.;
    }
    if (!face_is_wet()) {
      for (vector uf in ufl)
	uf.x[] = 0.;
    }
    else
#endif // TREE

    if (H + Hm > 0.) {
      H = (H + Hm)/2.;
      double cp = hydrostatic ? sqrt(G*H) : sqrt(G*Delta*tanh(H/Delta));
      cp /= CFL_H;
      for (vector uf in ufl) {
	double c = fm.x[]*cp + fabs(uf.x[]);
	if (c > 0.) {
	  double dt = cm[]*Delta/c;
	  if (dt < dtmax)
	    dtmax = dt;
	}
      }
    }
  }
  dt = dtnext (CFL*dtmax);
}

#ifdef F0
event coriolis (i++)
{
  foreach() {
    double cosomega = cos (F0()*dt);
    double sinomega = sin (F0()*dt);
    for (vector u in ul) {
      double ua = u.x[];
      u.x[] =   ua*cosomega + u.y[]*sinomega;
      u.y[] = - ua*sinomega + u.y[]*cosomega;
    }
  }
}
#endif

/**
Vertical diffusion (including viscosity) is added by this code. */

#include "diffusion.h"

/**
## Mesh adaptation and cleanup

Adaptation plugs itself here. The fields and lists allocated in
[`defaults()`](#defaults0) above must be freed at the end of the run. */

#if TREE
event adapt (i++,last);
#endif

event cleanup (i = end, last)
{
  for (int l = 0; l < nl; l++)
    free (tracers[l]);
  free (tracers), tracers = NULL;
  delete ((scalar *) ul), free (ul), ul = NULL;
  delete ((scalar *) ufl), free (ufl), ufl = NULL;
  delete ((scalar *) al), free (al), al = NULL;
  delete (hl), free (hl), hl = NULL;
  delete ({eta});
}
