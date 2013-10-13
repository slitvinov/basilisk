#include "hele-shaw.h"

#define MAXLEVEL 8

double mu1 = 2.5, mu2 = 1e-3, k = 0.1;

p[left]  = dirichlet(1e-3);
p[right] = dirichlet(0);
f[left]  = 1.;

void parameters()
{
  L0 = 7.5e-2;
  N = 1 << MAXLEVEL;
  TOLERANCE = 1e-3;
}

event init (i = 0) {
  foreach()
    f[] = (x < 1e-3)*(1. - 1e-3*(1. + noise()));
}

event logfile (i++)
{
  double flux = 0.;
  foreach_boundary (left, false)
    flux += u.x[]*Delta;
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g %g\n", 
	   i, t, mgp.i, flux, s.sum, s.min, s.max);
}

event movies (t += 0.2)
{
  static FILE * fp1 = popen ("ppm2mpeg > f.mpg", "w");
  output_ppm (f, fp1, N, min = 0, max = 1);

  static FILE * fp2 = popen ("ppm2mpeg > v.mpg", "w");
  output_ppm (u.y, fp2, N);

  static FILE * fp3 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp3, N, min = 0, max = MAXLEVEL, linear = false);
}

event output (t += 10; t <= 50)
{
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({p,f,u}, stdout, N);
}

#if QUADTREE
event adapt (i++) {
  double tolerance = 0.02;
  adapt_wavelet ({f}, &tolerance, MAXLEVEL, list = {f});
}
#endif

event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1,0])/2.;
    kmu.x[] = - k/(mu1 + clamp(ff,0,1)*(mu2 - mu1));
  }
}

int main() { run(); }
