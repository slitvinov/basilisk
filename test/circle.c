/* This is similar to gerris/test/poisson/circle */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "mg.h"
#include "adapt.h"

scalar a = new scalar, b = new scalar, res = new scalar, dp = new scalar;

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

void boundary_dirichlet (scalar v)
{
  /* Dirichlet condition on all boundaries */
  foreach_boundary (right)
    v[ghost] = 2.*solution(x + delta/2., y) - v[];
  foreach_boundary (left)
    v[ghost] = 2.*solution(x - delta/2., y) - v[];
  foreach_boundary (top)
    v[ghost] = 2.*solution(x, y + delta/2.) - v[];
  foreach_boundary (bottom)
    v[ghost] = 2.*solution(x, y - delta/2.) - v[];
  restriction (v, v);
  update_halo (-1, v, v);
}

void homogeneous_boundary (scalar v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (right, l)  v[ghost] = - v[];
  foreach_boundary_level (left, l)   v[ghost] = - v[];
  foreach_boundary_level (top, l)    v[ghost] = - v[];
  foreach_boundary_level (bottom, l) v[ghost] = - v[];
  /* we don't need to restrict because the solution is already defined
     on coarse levels */
  update_halo (l, v, v);
}

void relax (scalar a, scalar b, int l)
{
  foreach_level (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - delta*delta*b[])/4.;
}

void residual (scalar a, scalar b, scalar res)
{
#if 1
  /* conservative coarse/fine discretisation */
  vector g = new vector;
  foreach() {
    g.x[] = (a[] - a[-1,0])/delta;
    g.y[] = (a[] - a[0,-1])/delta;
  }
  foreach_boundary (right)
    g.x[1,0] = (a[1,0] - a[])/delta;
  foreach_boundary (top)
    g.y[0,1] = (a[0,1] - a[])/delta;
  restriction_u_v (g.x, g.y);
  update_halo_u_v (-1, g.x, g.y);
  foreach()
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/delta;
#else
  /* "naive" discretisation */
  foreach()
    res[] = b[] + 
    (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/(delta*delta);
#endif
}

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  QUADTREE_VARIABLES;
  VARIABLES;
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

void solve (int depth)
{
  int nrelax = 4;
  init_grid(1);

  while (refine_function (0, nvar - 1, refine_circle, &depth));
  flag_halo_cells();
  foreach()
    b[] = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  boundary_dirichlet (a);

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (a, res, dp,
	      relax, homogeneous_boundary,
	      nrelax, 0);
    boundary_dirichlet (a);
    residual (a, b, res);
    double max = 0.;
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach() {
    double e = a[] - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "max error %d %g\n", depth, max);

  free_grid();
}

int main (int argc, char ** argv)
{
  for (int depth = 7; depth <= 10; depth++)
    solve (depth);
}
