#include "poisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar alpha;
  double dt;
};

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = data;
  (const) face vector mu = p->mu;
  (const) scalar alpha = p->alpha;
  double dt = p->dt;
  vector u = {a[0], a[1]}, r = {b[0], b[1]};

  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (dt*alpha[]*(mu.x[1,0]*u.x[1,0] + mu.x[]*u.x[-1,0] +
			   mu.y[0,1]*u.x[0,1] + mu.y[]*u.x[0,-1]
#if IMPLICIT
			   + ((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
			      (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/2.
#endif
			   ) + r.x[]*sq(Delta))/
    (sq(Delta) + dt*alpha[]*(mu.x[1,0] + mu.x[] + mu.y[0,1] + mu.y[]));

#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = data;
  (const) face vector mu = p->mu;
  (const) scalar alpha = p->alpha;
  double dt = p->dt;
  vector u = {a[0], a[1]}, r = {b[0], b[1]}, res = {resl[0], resl[1]};
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  foreach_dimension() {
    face vector g[];
    scalar a = u.x;
    foreach_face()
      g.x[] = mu.x[]*(a[] - a[-1,0])/Delta;
    boundary_flux ({g});
    foreach (reduction(max:maxres)) {
      res.x[] = r.x[] - u.x[] + dt*alpha[]/Delta*
        (g.x[1,0] - g.x[] + g.y[0,1] - g.y[]
#if IMPLICIT
	 + ((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
	    (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/(2.*Delta)
#endif
	 );
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - u.x[] +
        dt*alpha[]*(mu.x[1,0]*u.x[1,0] + mu.x[]*u.x[-1,0] +
		    mu.y[0,1]*u.x[0,1] + mu.y[]*u.x[0,-1] -
		    (mu.x[1,0] + mu.x[] + mu.y[0,1] + mu.y[])*u.x[]
#if IMPLICIT
		    + ((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
		       (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/2.
#endif
		    )/sq(Delta);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
#endif
  return maxres;
}

mgstats viscosity (struct Viscosity p)
{
  vector u = p.u;
  (const) face vector mu = p.mu;
  (const) scalar alpha = p.alpha;
  double dt = p.dt;
  vector r[];
  foreach() {
    foreach_dimension()
      r.x[] = u.x[]
#if !IMPLICIT
      + dt*alpha[]*((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
		    (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/(2.*sq(Delta))
#endif
      ;
  }
  restriction ({mu,alpha});
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p);
}
