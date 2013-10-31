/**
# Bell-Collela-Glaz advection scheme

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89). Given a centered
scalar field `f`, a staggered vector field `u` and a timestep `dt`, it
fills the staggered vector field `flux` with the components of the
advection fluxes of `f`. */

void tracer_fluxes (const scalar f, 
		    const staggered vector u,
		    staggered vector flux,
		    double dt)
{

/**
We first compute the cell-centered gradient of `f` in a locally-allocated
vector field. */

  vector g[];
  gradients ({f}, {g});

/**
For each face, the flux is composed of two parts... */

  trash ({flux});
  foreach_face() {

/**
A normal component... */

    double un = dt*u.x[]/Delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*Delta/2.;

/**
and a tangential component... */

    double vn = u.y[i,0] + u.y[i,1];
    double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
    f2 -= dt*vn*fyy/(4.*Delta);
    flux.x[] = f2*u.x[];
  }

/**
Boundary conditions ensure the consistency of fluxes across
variable-resolution boundaries (on adaptive meshes). */

  boundary_normal ({flux});
}

/**
The function above is then used to do the time-integration itself,
using the generic tracer advection event. */

#include "tracer.h"
