#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

/* Dirichlet condition on all boundaries */
a[right]  = 2.*solution(x, y) - a[];
a[left]   = 2.*solution(x, y) - a[];
a[top]    = 2.*solution(x, y) - a[];
a[bottom] = 2.*solution(x, y) - a[];

int main (int argc, char ** argv)
{
  X0 = Y0 = -0.5;
  int depth = argc < 2 ? 9 : atoi(argv[1]), nrelax = 4;
  init_grid(1 << depth);

  foreach() {
    b[] = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
    a[] = 0.;
  }
  boundary ({a});

  #define NITER 13
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  residual (a, b, res);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (a, res, dp,
	      relax, homogeneous_boundary,
	      nrelax, 0);
    boundary ({a});
    residual (a, b, res);
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
