/* This is similar to gerris/test/poisson/circle */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "mg.h"
#include "adapt.h"

var a = new var, b = new var, res = new var, dp = new var;

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

void boundary (void * grid, var v)
{
  /* Dirichlet condition on all boundaries */
  foreach_boundary (grid, right)
    v(+1,0) = 2.*solution(x + delta/2., y) - v(0,0);
  foreach_boundary (grid, left)
    v(-1,0) = 2.*solution(x - delta/2., y) - v(0,0);
  foreach_boundary (grid, top)
    v(0,+1) = 2.*solution(x, y + delta/2.) - v(0,0);
  foreach_boundary (grid, bottom)
    v(0,-1) = 2.*solution(x, y - delta/2.) - v(0,0);
  restriction (grid, v);
  update_halo (grid, -1, v, v);
}

void homogeneous_boundary (void * grid, var v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (grid, right, l)   v(+1,0) = - v(0,0);
  foreach_boundary_level (grid, left, l)    v(-1,0) = - v(0,0);
  foreach_boundary_level (grid, top, l)     v(0,+1) = - v(0,0);
  foreach_boundary_level (grid, bottom, l)  v(0,-1) = - v(0,0);
  /* we don't need to restrict because the solution is already defined
     on coarse levels */
  update_halo (grid, l, v, v);
}

void relax (void * grid, var a, var b, int l)
{
  foreach_level (grid, l)
    a(0,0) = (a(1,0) + a(-1,0) + a(0,1) + a(0,-1) 
	      - delta*delta*b(0,0))/4.;
}

void residual (void * grid, var a, var b, var res)
{
  foreach (grid)
    res(0,0) = b(0,0) + 
    (4.*a(0,0) - a(1,0) - a(-1,0) - a(0,1) - a(0,-1))/(delta*delta);
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
  void * grid = init_grid(1);

  while (refine_function (grid, 1, 0, refine_circle, &depth));
  flag_halo_cells (grid);
  foreach(grid)
    b(0,0) = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  boundary (grid, a);

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (grid, a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (grid, a, res, dp,
	      relax, homogeneous_boundary,
	      nrelax, 0);
    boundary (grid, a);
    residual (grid, a, b, res);
    double max = 0.;
    foreach(grid)
      if (fabs(res(0,0)) > max)
	max = fabs(res(0,0));
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach(grid) {
    double e = a(0,0) - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a(0,0), b(0,0), res(0,0), e);
  }
  fprintf (stderr, "max error %d %g\n", depth, max);

  free_grid(grid);
}

int main (int argc, char ** argv)
{
  for (int depth = 7; depth <= 10; depth++)
    solve (depth);
}
