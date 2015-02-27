// the MPI parallel version of bump2D.c

#include "saint-venant.h"

#define LEVEL 7

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (i = 0)
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  if (pid() == 0)
    fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 2.5; t += 2.5/30) {
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, n = 256, min = 0, max = LEVEL);
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  if (pid() == 0)
    fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
