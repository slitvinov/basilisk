#include "grid/cartesian1D.h"
#include "saint-venant2.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (500);
  run();
}

void init()
{
  foreach() {
    double epsilon = 1e-2;
    w[] = 1. + epsilon*(-0.4 < x && x < -0.3);
    B[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
  }
}

event logfile (t = 0.7) {
  foreach()
    fprintf (stderr, "%g %g %.6f %g\n", x, w[], hu[], B[]);
}
