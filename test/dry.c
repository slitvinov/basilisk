#include "grid/cartesian1D.h"
#include "saint-venant1.h"

void parameters()
{
  N = 200;
}

void init()
{
  q.x[left] = - q.x[];
  h[right] = 0.;
  foreach() {
    zb[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
    h[] = 0.8 - zb[];
  }
}

event (t = {0.5, 0.75, 1, 3, 50}) {
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, h[], q.x[], zb[]);
  fprintf (stderr, "\n");
}

int main() { run(); }
