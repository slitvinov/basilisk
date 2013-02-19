#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include "mg.h"

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
}

void homogeneous_boundary (void * grid, var v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (grid, right, l)   v(+1,0) = - v(0,0);
  foreach_boundary_level (grid, left, l)    v(-1,0) = - v(0,0);
  foreach_boundary_level (grid, top, l)     v(0,+1) = - v(0,0);
  foreach_boundary_level (grid, bottom, l)  v(0,-1) = - v(0,0);
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

int main(int argc, char ** argv)
{
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  void * grid = init_grid(1 << depth);

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
    fprintf (stderr, "%d %g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }
  double max = 0;
  foreach(grid) {
    double e = a(0,0) - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a(0,0), b(0,0), res(0,0), e);
  }
  fprintf (stderr, "# max error %g\n", max);

  free_grid(grid);
}
