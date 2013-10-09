// Poisson equation with variable coefficients

#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];
vector c[];

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

/* Dirichlet condition on all boundaries */
a[right]  = dirichlet (solution(x, y));
a[left]   = dirichlet (solution(x, y));
a[top]    = dirichlet (solution(x, y));
a[bottom] = dirichlet (solution(x, y));

static void relax_variable (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l) {
    //    fprintf (stderr, "%d %g %g c: %g %g %g %g\n",
    //	     l, x, y, c.x[], c.x[1,0], c.y[], c.y[0,1]);
    a[] = (c.x[1,0]*a[1,0] + c.x[]*a[-1,0] + c.y[0,1]*a[0,1] + c.y[]*a[0,-1] 
	   - sq(Delta)*b[])/(c.x[1,0] + c.x[] + c.y[0,1] + c.y[]);
  }
}

static double residual_variable (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = (a[] - a[-1,0])/Delta;
  boundary_normal ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (c.x[]*g.x[] - c.x[1,0]*g.x[1,0] + 
		   c.y[]*g.y[] - c.y[0,1]*g.y[0,1])/Delta;
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + 
      ((c.x[1,0] + c.x[] + c.y[0,1] + c.y[])*a[]
       - c.x[1,0]*a[1,0] - c.x[]*a[-1,0] - c.y[0,1]*a[0,1] - c.y[]*a[0,-1])
      /sq(Delta);
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

int main (int argc, char ** argv)
{
  X0 = Y0 = -0.5;
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  init_grid(1 << depth);

  foreach() {
    b[] = -18.*sq(pi)*sin(3.*pi*x)*sin(3.*pi*y)*(x + y + 2.) +
      3.*pi*cos(3.*pi*x)*sin(3.*pi*y) +
      3.*pi*sin(3.*pi*x)*cos(3.*pi*y);
    a[] = 0.;
  }
  boundary ({a});

  trash ({c});
  foreach_face()
    c.x[] = x + y + 2.;
  foreach_fine_to_coarse()
    foreach_dimension()
      c.x[] = (fine(c.x,0,0) + fine(c.x,0,1))/2.;
  boundary_normal ({c});
  foreach_boundary_fine_to_coarse(right)
    c.x[] = (fine(c.x,0,0) + fine(c.x,0,1))/2.;
  foreach_boundary_fine_to_coarse(top)
    c.y[] = (fine(c.y,0,0) + fine(c.y,1,0))/2.;

  #define NITER 13
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual_variable (a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (a, res, dp, relax_variable, nrelax, 0);
    residual_variable (a, b, res);
    double max = 0.;
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "%d %.2g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, 
    	    maxres[i]);
  }
  double max = 0;
  foreach() {
    double e = a[] - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "# max error %g\n", max);

  free_grid();
}
