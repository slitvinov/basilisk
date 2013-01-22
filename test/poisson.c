#include <math.h>
#include <stdlib.h>
#include <time.h>

struct _Data {
  double a, b, res, dp;
};

#include GRID // works with all "multigrid" grids i.e. quadtree.h or multigrid.h
#include "utils.h"
#include "wavelet.h"

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
}

void boundary_level (void * grid, var v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (grid, right, l)   val(v,+1,0) = - val(v,0,0);
  foreach_boundary_level (grid, left, l)    val(v,-1,0) = - val(v,0,0);
  foreach_boundary_level (grid, top, l)     val(v,0,+1) = - val(v,0,0);
  foreach_boundary_level (grid, bottom, l)  val(v,0,-1) = - val(v,0,0);
}

void relax (void * grid, var a, var b, int l)
{
  foreach_level (grid, l)
    val(a,0,0) = (val(a,1,0) + val(a,-1,0) +
		  val(a,0,1) + val(a,0,-1) 
		  - delta*delta*val(b,0,0))/4.;
  boundary_level (grid, a, l);
}

void residual (void * grid, var a, var b, var res)
{
  foreach (grid)
    val(res,0,0) = val(b,0,0) + 
    (4.*val(a,0,0) - val(a,1,0) - val(a,-1,0) - val(a,0,1) - val(a,0,-1))/(delta*delta);
}

void cycle (void * grid, int depth, var a, var res, var dp, 
	    int nrelax, int minlevel)
{
  restriction (grid, res);
  for (int l = minlevel; l <= depth; l++) {
    if (l == minlevel) {
      foreach_level (grid, l)
	val(dp,0,0) = 0.;
    }
    else /* linear interpolation from coarser level (works better than bilinear) */
      foreach_level (grid, l)
	val(dp,0,0) = coarse(dp,0,0) + 
	  (childx*(coarse(dp,1,0) - coarse(dp,-1,0)) +
	   childy*(coarse(dp,0,1) - coarse(dp,0,-1)))/8.;
    boundary_level (grid, dp, l);
    for (int i = 0; i < nrelax; i++)
      relax (grid, dp, res, l);
  }
  foreach(grid)
    val(a,0,0) += val(dp,0,0);
  boundary (grid, a);
}

int main(int argc, char ** argv)
{
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  void * grid = init_grid(1 << depth);
  var a = var(a), b = var(b), res = var(res), dp = var(dp);

  foreach(grid)
    val(b,0,0) = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  boundary (grid, a);

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (grid, a, b, res);
  for (int i = 0; i < NITER; i++) {
    cycle (grid, depth, a, res, dp, nrelax, 0);
    residual (grid, a, b, res);
    double max = 0.;
    foreach(grid)
      if (fabs(val(res,0,0)) > max)
	max = fabs(val(res,0,0));
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "%d %g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach(grid) {
    double e = val(a,0,0) - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, val(a,0,0), val(b,0,0), val(res,0,0), e);
  }
  fprintf (stderr, "# max error %g\n", max);

  free_grid(grid);
}
