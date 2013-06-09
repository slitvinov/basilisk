#include "grid/cartesian1D.h"
#include "saint-venant.h"

void parameters()
{
  X0 = Y0 = -0.5;
  N = 200;
}

void init()
{
  h[right]   = 0.;
  eta[right] = zb[];
  u.x[right] = u.x[];
  foreach() {
    zb[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
    h[] = 0.8 - zb[];
  }
}

event logfile (t = {0.5, 0.75, 1, 3, 50}) {
  foreach() {
    fprintf (stderr, "%g %g %g %g\n", x, h[], u.x[], zb[]);
    assert (h[] > 0.);
  }
  fprintf (stderr, "\n");
}

int main() { run(); }
