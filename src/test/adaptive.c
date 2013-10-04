// similar to advection.c but with adaptivity

#include "advection.h"

void parameters()
{
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
}

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

void init()
{
  foreach()
    f[] = bump(x,y);
}

event output (t++) {
  static int nf = 0;
  printf ("file: f-%d\n", nf);
  output_field ({f}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, stdout, N);
}

event velocity (i++) {
  adapt_wavelet ({f}, (double[]){1e-2}, 8, 5, list = {f});

  trash ({u});
  foreach_face(x)
    u.x[] = 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*cos((y + 0.5)*pi);
  foreach_face(y)
    u.y[] = - 1.5*sin(2.*pi*t/5.)*cos((x + 0.5)*pi)*sin((y + 0.5)*pi);
  boundary_mac ({u});
}

event logfile (t = {0,5}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
  static double sum = 0.;
  if (t > 0.)
    assert (fabs(s.sum - sum) < 1e-10); // tracer conservation
  sum = s.sum;
}

event field (t = 5) {
  scalar e[];
  foreach()
    e[] = f[] - bump(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  
  if (N == 256) {
    printf ("file: error\n");
    output_field ({e}, stdout, N, linear = false);
  }
}

int main() {
  for (N = 256; N <= 256; N *= 2)
    run();
}
