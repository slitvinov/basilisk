// fixme: utils.h should be in stream.h
// #include "utils.h"
#include "navier-stokes/stream.h"

#define MAXLEVEL 8

void parameters()
{
  X0 = Y0 = -0.5;
  N = 1 << MAXLEVEL;
}

void init()
{
  double dd = 0.1;
  foreach()
    omega[] = (exp(-(sq(x - dd) + sq(y))/(dd/10.)) +
	       exp(-(sq(x + dd) + sq(y))/(dd/10.)));
}

event logfile (t = {0,30}) {
  stats s = statsf (omega);
  fprintf (stderr, "%g %d %g %g %d\n", t, i, dt, s.sum, mgpsi.i);
}

event output (t += 5) {
  static int nf = 0;
  printf ("file: psi-%d\n", nf);
  output_field ({psi, omega}, stdout, N, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({l}, stdout, N);
  nf++;
}

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({omega}, (double[]){1e-2}, MAXLEVEL, list = {omega, psi});
}
#endif

int main() { run(); }
