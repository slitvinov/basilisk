#include "saint-venant1.h"

#define LEVEL 7

void parameters()
{
  X0 = Y0 = -0.5;
  N = 1 << LEVEL;
  G = 9.81;
  L0 = 1000.;
}

#define zb(x,y) ((cos(pi*x/L0)*cos(pi*y/L0) + \
		  cos(3.*pi*x/L0)*cos(3.*pi*y/L0)) - 2.*x/1000.)

void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = zb(x,y);
}

void init()
{
  zb.refine = refine_zb; // updates terrain
  foreach() {
    zb[] = zb(x,y);
    h[] = 0.1;
  }
}

int event (i++) {
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
  boundary ((scalar *){u});
}

int event (i += 10) {
  stats s = statsf (h);
  norm n = normf (u.x);
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g\n", 
	   t, i, s.min, s.max, s.sum, dt, n.rms, n.max);
  assert (s.min > 0.);
}

int event (t <= 1200.; t += 1200./8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h, zb, u}, N, stdout, true);

  scalar l = new scalar;
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({h, zb, l}, N, stdout, false);
  delete ({l});

  nf++;
}

// int event (t += 100)
//  output_matrix (h, N, stdout, true);

int event (i++) {
  scalar eta = new scalar, w = new scalar;
  foreach() eta[] = h[] + zb[];
  boundary ({eta});
  wavelet (h, w);

  double cmax = 1e-2;
  scalar * list = {h, zb, u, dh, dq};
  int nf = refine_wavelet (w, cmax, LEVEL, list);
  int nc = coarsen_wavelet (w, cmax/4., 4, list);
  if (nf || nc)
    boundary (list);
  delete ({eta, w});

  //  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int main() { run(); }
