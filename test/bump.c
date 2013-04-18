#include "grid/cartesian1D.h"
#include "saint-venant.h"

void parameters()
{
  N = 500;
}

void init()
{
  foreach() {
    double epsilon = 1e-2;
    w[] = 1. + epsilon*(-0.4 < x && x < -0.3);
    B[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
  }
}

int event (t = 0.7) {
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, w[], hu[], B[]);
}

int main() { run(); }
