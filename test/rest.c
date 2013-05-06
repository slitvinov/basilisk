// lake at rest with variable resolution
#include "saint-venant1.h"

void parameters()
{
  N = 16;
}

int refine (Point point, void * data)
{
  return x < 0.1;
}

void init()
{
  refine_function (refine, NULL, scalars (h, zb, q, dh, dq));

  foreach() {
    zb[] = 0.2*exp(-200*x*x);
    h[] = 1. - zb[];
  }
}

int event (i = 1)
{
  norm n = normf (q.x);
  fprintf (stderr, "%g %g %g\n", n.avg, n.rms, n.max);
  foreach ()
    printf ("%g %g %g %g %g\n", x, y, h[], zb[], q.x[]);
}

int main() { run(); }
