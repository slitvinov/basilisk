#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (200);
  run();
}

h[right]   = 0.;
eta[right] = zb[];
u.n[right] = u.x[];

event init (i = 0)
{
  foreach() {
    zb[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
    h[] = 0.8 - zb[];
  }
}

event logfile (t = {0.5, 0.75, 1, 3, 50}) {
  foreach() {
    fprintf (stderr, "%g %g %.6f %g\n", x, h[], u.x[], zb[]);
    assert (h[] > 0.);
  }
  fprintf (stderr, "\n");
}
