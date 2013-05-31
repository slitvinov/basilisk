// Small amplitude solitary wave interacting with a parabolic hump
#include "saint-venant1.h"

#define MAXLEVEL 9
#define MINLEVEL 7

void parameters()
{
  N = 1 << MAXLEVEL;
  X0 = -0.5;
  Y0 = -1.;
  L0 = 2.;
}

u.x[left]   = u.x[];
u.x[right]  = u.x[];
u.y[top]    = u.y[];
u.y[bottom] = u.y[];

void init()
{
  foreach() {
    x += 0.5;
    zb[] = 0.8*exp(-5.*sq(x - 0.9) - 50.*y*y);
    h[] = (0.05 < x && x < 0.15 ? 1.01 : 1) - zb[];
  }
}

event outputfile (t = {0.6, 0.9, 1.2, 1.5, 1.8})
{
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, N, stdout, true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, N, stdout, false);
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event adapt (i++) {
  scalar w[];
  wavelet (eta, w);

  double cmax = 1e-4;
  int nf = refine_wavelet (w, cmax, MAXLEVEL, all);
  int nc = coarsen_wavelet (w, cmax/4., MINLEVEL, all);
  if (nf || nc)
    boundary (all);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int main() { run(); }
