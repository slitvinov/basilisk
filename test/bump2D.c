#include "saint-venant.h"

#define LEVEL 7

void parameters()
{
  X0 = Y0 = -0.5;
  N = 1 << LEVEL;
}

void init()
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 2.5; t += 2.5/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, stdout, N);

  /* check symmetry */
  foreach() {
    double h0 = h[];
    point = locate (-x, -y);
    //    printf ("%g %g %g %g %g\n", x, y, h0, h[], h0 - h[]);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (-x, y);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (x, -y);
    assert (fabs(h0 - h[]) < 1e-12);
  }
}

event adapt (i++) {
  scalar w[];
  wavelet (h, w);

  double cmax = 1e-3;
  int nf = refine_wavelet (w, cmax, LEVEL, all);
  int nc = coarsen_wavelet (w, cmax/4., 0, all);
  if (nf || nc)
    boundary (all);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int main() { run(); }
