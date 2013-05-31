#include "grid/cartesian1D.h"
#include "saint-venant1.h"

#define LEVEL 10

void parameters()
{
  X0 = -0.5;
  N = 1 << LEVEL;
}

void init()
{
  foreach() {
    h[] = exp(-200.*(x*x));
    x -= -0.25;
    zb[] = exp(-200.*(x*x));
    h[] = max(h[] - zb[], 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f %g\n", t, i, s.min, s.max, s.sum, dt);
  assert (s.min >= 0.);
}

event outputfile (t <= 0.6; t += 0.6/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf++);
  foreach()
    printf ("%g %g %g %g\n", x, h[], zb[], u.x[]);
}

int main() { run(); }
