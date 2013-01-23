/* This is similar to gerris/test/poisson/circle */

#include <math.h>
#include <stdlib.h>
#include <time.h>

struct _Data {
  double a, b, res, dp;
};

#include "grid/quadtree.h"
#include "utils.h"
#include "mg.h"
#include "adapt.h"

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

void boundary (void * grid, var v)
{
  /* Dirichlet condition on all boundaries */
  foreach_boundary (grid, right)
    val(v,+1,0) = 2.*solution(x + delta/2., y) - val(v,0,0);
  foreach_boundary (grid, left)
    val(v,-1,0) = 2.*solution(x - delta/2., y) - val(v,0,0);
  foreach_boundary (grid, top)
    val(v,0,+1) = 2.*solution(x, y + delta/2.) - val(v,0,0);
  foreach_boundary (grid, bottom)
    val(v,0,-1) = 2.*solution(x, y - delta/2.) - val(v,0,0);
  restriction (grid, v);
  update_halos (grid, -1, v, v);
}

void homogeneous_boundary (void * grid, var v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (grid, right, l)   val(v,+1,0) = - val(v,0,0);
  foreach_boundary_level (grid, left, l)    val(v,-1,0) = - val(v,0,0);
  foreach_boundary_level (grid, top, l)     val(v,0,+1) = - val(v,0,0);
  foreach_boundary_level (grid, bottom, l)  val(v,0,-1) = - val(v,0,0);
  /* we don't need to restrict because the solution is already defined
     on coarse levels */
  update_halos (grid, l, v, v);
}

void relax (void * grid, var a, var b, int l)
{
  foreach_level (grid, l)
    val(a,0,0) = (val(a,1,0) + val(a,-1,0) +
		  val(a,0,1) + val(a,0,-1) 
		  - delta*delta*val(b,0,0))/4.;
}

void residual (void * grid, var a, var b, var res)
{
  foreach (grid)
    val(res,0,0) = val(b,0,0) + 
    (4.*val(a,0,0) - val(a,1,0) - val(a,-1,0) - val(a,0,1) - val(a,0,-1))/(delta*delta);
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
  var a = var(a), b = var(b), res = var(res), dp = var(dp);
  int nrelax = 4;
  void * grid = init_grid(1);

  while (refine_function (grid, 1, 0, refine_circle, &depth));
  flag_halo_cells (grid);
  foreach(grid)
    val(b,0,0) = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  boundary (grid, a);

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (grid, a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (grid, depth, a, res, dp,
	      relax, homogeneous_boundary,
	      nrelax, 0);
    boundary (grid, a);
    residual (grid, a, b, res);
    double max = 0.;
    foreach(grid)
      if (fabs(val(res,0,0)) > max)
	max = fabs(val(res,0,0));
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach(grid) {
    double e = val(a,0,0) - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, val(a,0,0), val(b,0,0), val(res,0,0), e);
  }
  fprintf (stderr, "max error %d %g\n", depth, max);

  free_grid(grid);
}

int main (int argc, char ** argv)
{
  for (int depth = 7; depth <= 10; depth++)
    solve (depth);
}
