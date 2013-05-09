#include "grid/cartesian1D.h"
#include "saint-venant.h"

void parameters()
{
  X0 = Y0 = -0.5;
  N = 100;
}

void init()
{
  foreach() {
    w[] = 1.;
    B[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
    hu[] = 0.3*(w[] - B[]);
  }
}

int event (t = 1.8) {
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, w[], hu[], B[]);
}

int main() { run(); }
