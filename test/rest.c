// lake at rest with variable resolution
#include "saint-venant1.h"

void parameters()
{
  X0 = Y0 = -0.5;
  N = 16;
}

int refine (Point point, void * data)
{
  return x < 0.1 && y < 0.1;
}

void init()
{
  refine_function (refine, NULL, all);

  foreach() {
    zb[] = 0.2*exp(-200*(x*x + y*y));
    h[] = 1. - zb[];
  }
}

int event (i = 1)
{
  norm n = normf (u.x);
  fprintf (stderr, "# %g %g %g\n", n.avg, n.rms, n.max);
  foreach ()
    printf ("%g %g %g %g %g\n", x, y, h[], zb[], u.x[]);
}

int main() { run(); }
