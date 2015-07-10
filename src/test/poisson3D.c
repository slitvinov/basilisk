// similar to poisson.c but in 3D

#include "grid/multigrid3D.h"
#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];

double solution (double x, double y, double z)
{
  return sin(3.*pi*x)*sin(3.*pi*y)*sin(3.*pi*z);
}

int main (int argc, char ** argv)
{
  /* Dirichlet condition on all boundaries */
  foreach_dimension() {
    a[right] = dirichlet (solution(x, y, z));
    a[left]  = dirichlet (solution(x, y, z));
  }
  /* homogeneous conditions for dp */
  foreach_dimension() {
    dp[right] = dirichlet(0);
    dp[left]  = dirichlet(0);
  }
    
  origin (-0.5, -0.5, -0.5);
  int depth = argc < 2 ? 6 : atoi(argv[1]), nrelax = 4;
  init_grid(1 << depth);

  foreach() {
    b[] = -27.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y)*sin(3.*pi*z);
    a[] = 0.;
  }
  boundary ({a});

  #define NITER 13
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  const vector alpha[] = {1.,1.,1.};
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = alpha; p.lambda = lambda;
  scalar * lres = {res};
  residual ({a}, {b}, lres, &p);
  for (int i = 0; i < NITER; i++) {
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0);
    maxres[i] = residual ({a}, {b}, lres, &p);
    iter[i] = clock();
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "%d %.2g\n", i, maxres[i]);
    printf ("%d %g %g\n", i, (iter[i] - start)/(double)CLOCKS_PER_SEC, 
    	    maxres[i]);
  }
  double max = 0;
  foreach() {
    double e = a[] - solution(x, y, z);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "# max error %g\n", max);

  free_grid();
}
