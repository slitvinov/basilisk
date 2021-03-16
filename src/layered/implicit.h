#include "poisson.h"

mgstats mgH;
double theta_H = 0.5;

/**
## Setup

The $w_k$ and $\phi_k$ scalar fields are allocated and the $w_k$ are
added to the list of advected tracers. */

event defaults0 (i = 0)
{
  hydrostatic = true;
  if (CFL_H == 1e40)
    CFL_H = HUGE;
  mgH.nrelax = 4;
}

trace
static void relax_hydro (scalar * ql, scalar * rhsl, int lev, void * data)
{
  scalar eta = ql[0], rhs_eta = rhsl[0];
  face vector alpha = *((vector *)data);
  foreach_level_or_leaf (lev) {
    double d = - cm[]*sq(Delta);
    double n = d*rhs_eta[];
    foreach_dimension() {
      n += alpha.x[1]*eta[1] + alpha.x[]*eta[-1];
      d += alpha.x[1] + alpha.x[];
    }
    eta[] = n/d;
  }
}

trace
static double residual_hydro (scalar * ql, scalar * rhsl,
			      scalar * resl, void * data)
{
  scalar eta = ql[0], rhs_eta = rhsl[0], res_eta = resl[0];
  face vector alpha = *((vector *)data);
  double maxres = 0.;
#if 1 // TREE
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(eta[] - eta[-1])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    /* conservative coarse/fine discretisation (2nd order) */
    res_eta[] = rhs_eta[] - eta[];
    foreach_dimension()
      res_eta[] -= (g.x[1] - g.x[])/(Delta*cm[]);
#endif // TREE

    if (fabs(res_eta[]) > maxres)
      maxres = fabs(res_eta[]);
  }
  boundary (resl);
  return maxres;
}

scalar res_eta[];

scalar rhs_eta;
face vector alpha_eta;

event acceleration (i++)
{
  if (theta_H < 1.)
    advect ((1. - theta_H)*dt);
    
  face vector su[];
  alpha_eta = new face vector;
  double C = - G*sq(theta_H*dt);
  foreach_face() {
    double ax = - theta_H*gmetric(0)*G*(eta[] - eta[-1])/Delta;
    su.x[] = alpha_eta.x[] = 0.;
    foreach_layer() {
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] = (1. - theta_H)*(hu.x[] + dt*hf.x[]*ax) + theta_H*hf.x[]*uf;
      ha.x[] -= hf.x[]*ax;
      su.x[] += hu.x[];
      alpha_eta.x[] += hf.x[];
    }
    alpha_eta.x[] *= C*gmetric(0);
  }
  boundary ((scalar *){alpha_eta, hu, hf, su, ha});

  rhs_eta = new scalar;
  foreach () {
    rhs_eta[] = eta[];
    foreach_dimension()
      rhs_eta[] -= dt*(su.x[1] - su.x[])/(Delta*cm[]);
  }

  // fixme: check which are really necessary
  restriction ({h, hf, alpha_eta, zb});
  
#if TREE
  eta.restriction = restriction_average;
#endif
}

event pressure (i++)
{
  mgH = mg_solve ({eta}, {rhs_eta}, residual_hydro, relax_hydro, &alpha_eta,
		  res = {res_eta},
		  nrelax = 4, minlevel = 1,
		  tolerance = TOLERANCE*(G*sq(dt)));
  delete ({rhs_eta, alpha_eta});
  
#if TREE
  eta.restriction = restriction_eta;
#endif

  foreach_face() {
    double ax = - theta_H*gmetric(0)*G*(eta[] - eta[-1])/Delta;
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] = theta_H*(hf.x[]*uf + dt*ha.x[]) - dt*ha.x[];
    }
  }
  boundary ((scalar *){ha, hu});
}
