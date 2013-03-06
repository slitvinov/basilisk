#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include "wavelet.h"
#include "mg.h"

scalar a = new scalar, b = new scalar, res = new scalar, dp = new scalar;

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

void boundary (scalar v)
{
  /* Dirichlet condition on all boundaries */
  foreach_boundary (right)
    v[+1,0] = 2.*solution(x + delta/2., y) - v[];
  foreach_boundary (left)
    v[-1,0] = 2.*solution(x - delta/2., y) - v[];
  foreach_boundary (top)
    v[0,+1] = 2.*solution(x, y + delta/2.) - v[];
  foreach_boundary (bottom)
    v[0,-1] = 2.*solution(x, y - delta/2.) - v[];
}

void homogeneous_boundary (scalar v, int l)
{
  /* Homogeneous Dirichlet condition on all boundaries */
  foreach_boundary_level (right, l)   v[+1,0] = - v[];
  foreach_boundary_level (left, l)    v[-1,0] = - v[];
  foreach_boundary_level (top, l)     v[0,+1] = - v[];
  foreach_boundary_level (bottom, l)  v[0,-1] = - v[];
}

void relax (scalar a, scalar b, int l)
{
  foreach_level (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - delta*delta*b[])/4.;
}

void residual (scalar a, scalar b, scalar res)
{
  foreach()
    res[] = b[] + 
    (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/(delta*delta);
}

int main(int argc, char ** argv)
{
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  init_grid(1 << depth);

  foreach()
    b[] = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  boundary (a);

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (a, res, dp,
	      relax, homogeneous_boundary,
	      nrelax, 0);
    boundary (a);
    residual (a, b, res);
    double max = 0.;
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "%d %g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
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
