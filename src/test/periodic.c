#include "poisson.h"
#include "utils.h"

int main()
{
  scalar a[], b[], e[];

  a[right] = periodic();
  a[top]   = periodic();

  for (int n = 8; n <= 128; n *= 2) {
    init_grid (n);
    foreach() {
      a[] = 0.;
      b[] = 8.*sq(pi)*sin(2.*pi*x)*cos(2.*pi*y);
    }
    boundary ({a,b});
    TOLERANCE = 1e-10;
    poisson (a, b);
    foreach()
      e[] = a[] + sin(2.*pi*x)*cos(2.*pi*y);
    stats s = statsf(e);
    foreach()
      e[] -= s.sum/s.area;
    fprintf (stderr, "%d %g\n", n, statsf(e).max);
  }
  foreach()
    printf ("%g %g %g %g\n", x, y, a[], e[]);
}
