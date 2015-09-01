#include "poisson.h"
#include "utils.h"

int main()
{
  scalar a[], b[], e[];

  foreach_dimension() {
    a[right] = periodic();
  }
  
  for (int n = 8; n <= 128; n *= 2) {
    init_grid (n);
    foreach() {
      a[] = 0.;
      b[] = 4.*dimension*sq(pi)*sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z);
    }
    boundary ({a,b});
    TOLERANCE = 1e-10;
    poisson (a, b);
    foreach()
      e[] = a[] + sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z);
    stats s = statsf(e);
    foreach()
      e[] -= s.sum/s.volume;
    fprintf (stderr, "%d %g\n", n, statsf(e).max);
  }
  foreach()
    printf ("%g %g %g %g %g\n", x, y, z, a[], e[]);
}
