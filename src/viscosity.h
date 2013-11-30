#include "poisson.h"

struct Viscosity {
  // mandatory
  vector u;
  double dt;
  face vector mu;
  vector r;
  // optional
  scalar alpha; // default 1
};

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = data;
  (const) face vector mu = p->mu;
  (const) scalar alpha = p->alpha;
  vector u = {a[0], a[1]}, r = {b[0], b[1]};

  foreach_level_or_leaf (l) {
    Delta *= Delta;
    foreach_dimension()
      u.x[] = (alpha[]*(mu.x[1,0]*u.x[1,0] + mu.x[]*u.x[-1,0] +
			mu.y[0,1]*u.x[0,1] + mu.y[]*u.x[0,-1] +
			((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
			 (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/2.) +
	       r.x[]*Delta)/
      (Delta + alpha[]*(mu.x[1,0] + mu.x[] + mu.y[0,1] + mu.y[]));      
  }
}

static double residual_viscosity (scalar * a, scalar * b, scalar ** presl, 
				  void * data)
{
  scalar * resl = *presl;

  if (!resl) {
    for (scalar s in a) {
      scalar res = new scalar;
      resl = list_append (resl, res);
    }
    *presl = resl;
  }

  struct Viscosity * p = data;
  (const) face vector mu = p->mu;
  (const) scalar alpha = p->alpha;
  vector u = {a[0], a[1]}, r = {b[0], b[1]}, res = {resl[0], resl[1]};
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  for (a,b,res in al,bl,resl) {
    vector g[];
    foreach_face()
      g.x[] = mu.x[]*(a[] - a[-1,0])/Delta;
    boundary_normal ({g});
    foreach (reduction(max:maxres)) {
      res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta
	- alpha[]*a[];
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    Delta *= Delta;
    foreach_dimension()
      res.x[] = r.x[] - u.x[] +
        alpha[]*(mu.x[1,0]*u.x[1,0] + mu.x[]*u.x[-1,0] +
		 mu.y[0,1]*u.x[0,1] + mu.y[]*u.x[0,-1] -
		 (mu.x[1,0] + mu.x[] + mu.y[0,1] + mu.y[])*u.x[] +
		 ((mu.x[1,0] - mu.x[])*(u.x[1,0] - u.x[-1,0]) + 
		  (mu.y[0,1] - mu.y[])*(u.y[1,0] - u.y[-1,0]))/2.)/Delta;
    if (fabs (res.x[]) > maxres)
      maxres = fabs (res.x[]);
  }
#endif
  return maxres;
}

mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r = p.r;

  foreach()
    foreach_dimension()
      r.x[] = u.x[] + dt*r.x[];

  if (p.alpha) {
    scalar alpha = p.alpha;
    foreach()
      alpha[] *= dt;
    restriction ({alpha});
  }
  else {
    const scalar alpha[] = dt;
    p.alpha = alpha;
  }

  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p);
}
