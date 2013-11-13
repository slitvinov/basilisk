/* This is similar to gerris/test/poisson/circle */

#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

/* Dirichlet condition on all boundaries */
a[right]  = dirichlet (solution(x, y));
a[left]   = dirichlet (solution(x, y));
a[top]    = dirichlet (solution(x, y));
a[bottom] = dirichlet (solution(x, y));

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

void solve (int depth)
{
  X0 = Y0 = -0.5;
  int nrelax = 4;
  init_grid(1);

  while (refine_function (refine_circle, &depth, NULL));
  foreach() {
    a[] = 0.;
    b[] = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  }
  boundary ({a});

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  const vector alpha[] = {1.,1.};
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = alpha; p.lambda = lambda;
  residual (a, b, res, &p);
  for (int i = 0; i < NITER; i++) {
    mg_cycle (a, res, dp, relax, &p, nrelax, 0);
    residual (a, b, res, &p);
    double max = 0.;
    foreach()
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %.1g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, 
	    (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
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
