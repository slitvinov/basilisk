#include "saint-venant1.h"

#define LEVEL 7

void parameters()
{
  X0 = Y0 = -0.5;
  N = 1 << LEVEL;
}

u.x[right]  = u.x[];
u.x[left]   = u.x[];
u.y[top]    = u.y[];
u.y[bottom] = u.y[];

void init()
{
  foreach() {
    h[] = exp(-200.*(x*x + y*y));
    x -= -0.25; y -= -0.25;
    zb[] = 0.5*exp(-200.*(x*x + y*y));
    h[] = max(h[] - zb[], 0.);
  }
}

int event (i++) {
  foreach() {
    x -= -0.25; y -= -0.25;
    zb[] = 0.5*exp(-200.*(x*x + y*y));
  }
  boundary ({zb});
}

int event (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);
  assert (s.min >= 0.);
}

int event (t <= 1.2; t += 1.2/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h}, N, stdout, true);

  scalar l = new scalar;
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, N, stdout, false);
}

int event (i++) {
  scalar w = new scalar;
  wavelet (h, w);

  double cmax = 1e-3;
  scalar * list = {h, zb, u, dh, dq};
  int nf = refine_wavelet (w, cmax, LEVEL, list);
  int nc = coarsen_wavelet (w, cmax/4., 0, list);
  if (nf || nc)
    boundary (list);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int main() { run(); }
