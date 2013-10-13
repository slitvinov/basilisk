/* Multigrid Poisson solvers */

static void boundary_homogeneous (scalar a, scalar da, int l)
{
  /* Apply homogeneous boundary conditions of a to da on all boundaries */
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true)
      da[ghost] = a.boundary_homogeneous[b] (point, da);
#if QUADTREE
  /* we don't need to restrict because the solution is already defined
     on coarse levels */
  halo_prolongation (l, {da});
#endif
}

void mg_cycle (scalar a, scalar res, scalar da,
	       void (* relax) (scalar da, scalar res, int depth),
	       int nrelax, int minlevel)
{
  restriction ({res});
  for (int l = minlevel; l <= depth(); l++) {
    if (l == minlevel) {
      foreach_level_or_leaf (l)
	da[] = 0.;
    }
    else 
      /* bilinear interpolation from coarser level */
      foreach_level_or_leaf (l)
	da[] = (9.*coarse(da,0,0) + 
		3.*(coarse(da,child.x,0) + coarse(da,0,child.y)) + 
		coarse(da,child.x,child.y))/16.;
    boundary_homogeneous (a, da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l);
      boundary_homogeneous (a, da, l);
    }
  }
  foreach()
    a[] += da[];
  boundary ({a});
}

// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on residual
double TOLERANCE = 1e-3;

typedef struct {
  int i;
  double maxres, sum;
} mgstats;

mgstats mg_solve (scalar a, scalar b, 
		  double (* residual) (scalar, scalar, scalar),
		  void (* relax) (scalar da, scalar res, int depth))
{
  scalar res[], da[];
  mgstats s = {0, 0., 0.};
  foreach () // fixme: need reduction(+:s.sum)
    s.sum += b[];
  s.maxres = residual (a, b, res);
  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.maxres > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, 4, 0);
    s.maxres = residual (a, b, res);
  }
  if (s.i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  sum: %g\n", 
	     NITERMAX, s.sum);
  return s;
}

// Poisson equation with constant coefficients

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

mgstats poisson (scalar a, scalar b)
{
  return mg_solve (a, b, residual, relax);
}

// Poisson equation with variable coefficients

vector alpha;

static void relax_variable (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (alpha.x[1,0]*a[1,0] + alpha.x[]*a[-1,0] + 
	   alpha.y[0,1]*a[0,1] + alpha.y[]*a[0,-1] 
	   - sq(Delta)*b[])
    /(alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[]);
}

static double residual_variable (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1,0])/Delta;
  boundary_normal ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta;
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + 
      ((alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[])*a[]
       - alpha.x[1,0]*a[1,0] - alpha.x[]*a[-1,0] 
       - alpha.y[0,1]*a[0,1] - alpha.y[]*a[0,-1])
      /sq(Delta);
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

mgstats poisson_variable (scalar a, scalar b, vector c)
{
  alpha = c;
  restriction_staggered ({alpha});
  return mg_solve (a, b, residual_variable, relax_variable);
}
