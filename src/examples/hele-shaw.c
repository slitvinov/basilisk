#include "grid/multigrid.h"
#include "hele-shaw.h"

double mu1 = 2.5, mu2 = 1e-3, k = 0.1;

p[left]  = dirichlet(1e-3);
p[right] = dirichlet(0);
f[left]  = 1.;

void parameters()
{
  L0 = 7.5e-2;
  N = 256;
  TOLERANCE = 1e-6;
}

void init() {
  foreach()
    f[] = (x < 1e-3)*(1. + 0.001*(1. + noise()));
}

event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1,0])/2.;
    kmu.x[] = - k/(mu1 + clamp(ff,0,1)*(mu2 - mu1));
  }
}

event logfile (i++)
{
  double flux = 0.;
  foreach_boundary (left, false)
    flux += u.x[]*Delta;
  fprintf (stderr, "%d %g %d %g\n", i, t, mgp.i, flux);
}

event movies (t += 0.2)
{
  static FILE * fp1 = NULL;
  if (!fp1) fp1 = popen ("ppm2mpeg > f.mpg", "w");
  output_ppm (f, fp1, min = 0, max = 1, n = N);

  static FILE * fp2 = NULL;
  if (!fp2) fp2 = popen ("ppm2mpeg > u.mpg", "w");
  output_ppm (u.x, fp2, n = N);

  static FILE * fp3 = NULL;
  if (!fp3) fp3 = popen ("ppm2mpeg > v.mpg", "w");
  output_ppm (u.y, fp3, n = N);
}

event output (t += 10; t <= 60)
{
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({p,f,u}, stdout, N);
}

int main() { run(); }
