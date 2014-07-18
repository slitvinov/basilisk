// lake at rest with emerged island (same as rest1.c but for Green-Naghdi)
#include "green-naghdi.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (32);
  run();
}

int refine (Point point, void * data)
{
  return x < 0.1 && y < 0.1;
}

event init (i = 0)
{
  refine_function (refine, NULL, all);

  foreach() {
    zb[] = 0.3*exp(-100*(x*x + y*y));
    h[] = max (0.1 - zb[], 0.);
  }
}

event logfile (i = 1)
{
  norm n = normf (u.x);
  fprintf (stderr, "# %.10f %.10f %.10f\n", n.avg, n.rms, n.max);

  FILE * fp = popen ("gfsview2D -s", "w");
  output_gfs (fp);
}
