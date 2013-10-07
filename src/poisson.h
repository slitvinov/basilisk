/* Multigrid Poisson solver */

void mg_cycle (scalar a, scalar res, scalar dp,
	       void (* relax)    (scalar dp, scalar res, int depth),
	       void (* boundary) (scalar * dp, int depth),
	       int nrelax, int minlevel)
{
  restriction ({res});
  for (int l = minlevel; l <= depth(); l++) {
    if (l == minlevel) {
      foreach_level_or_leaf (l)
	dp[] = 0.;
    }
    else 
      /* bilinear interpolation from coarser level */
      foreach_level_or_leaf (l)
	dp[] = (9.*coarse(dp,0,0) + 
		3.*(coarse(dp,child.x,0) + coarse(dp,0,child.y)) + 
		coarse(dp,child.x,child.y))/16.;
    boundary ({dp}, l);
    for (int i = 0; i < nrelax; i++) {
      relax (dp, res, l);
      boundary ({dp}, l);
    }
  }
  foreach()
    a[] += dp[];
}

// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on residual
double TOLERANCE = 1e-3;

static void homogeneous_boundary (scalar * v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  scalar p = *v;
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true)
      p[ghost] = - p[];
#if QUADTREE
  /* we don't need to restrict because the solution is already defined
     on coarse levels */
  halo_prolongation (l, {p});
#endif
}

static void relax (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - sq(Delta)*b[])/4.;
}

static double residual (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = (a[] - a[-1,0])/Delta;
  boundary_normal ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta;
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/sq(Delta);
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

typedef struct {
  int i;
  double maxres, sum;
} mgstats;

mgstats poisson (scalar a, scalar b)
{
  scalar res[], da[];
  mgstats s;
  foreach (reduction(+:sum))
    s.sum += b[];
  s.maxres = residual (a, b, res);
  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.maxres > TOLERANCE); s.i++) {
    mg_cycle (a, res, da,
	      relax, homogeneous_boundary,
	      4, 0);
    boundary ({a});
    s.maxres = residual (a, b, res);
  }
  if (s.i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  sum: %g\n", 
	     NITERMAX, s.sum);
  return s;
}
