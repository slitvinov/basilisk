#include "saint-venant1.h"

#define LEVEL 7

void parameters()
{
  N = 1 << LEVEL;
}

void init()
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

scalar w = new scalar;

int event (i++) {

  wavelet (h, w);

  double cmax = 1e-3;
  scalar * list = scalars (h, zb, q, dh, dq);
  int nf = refine_wavelet (w, cmax, LEVEL, list);
  int nc = coarsen_wavelet (w, cmax/4., 0, list);
  if (nf || nc)
    boundary (list);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int event (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

int event (t <= 2.5; t += 2.5/8) {
  scalar eta = new scalar;
  foreach()
    eta[] = h[] > 1e-3 ? h[] + zb[] : 0.;
  boundary (eta);
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field (eta, N, stdout, true);

  scalar l = new scalar;
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field (l, N, stdout, false);

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

int main() { run(); }
