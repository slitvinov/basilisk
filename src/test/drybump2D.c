#include "saint-venant.h"

#define LEVEL 7

int main()
{
  X0 = Y0 = -0.5;
  N = 1 << LEVEL;
  run();
}

u.x[right]  = u.x[];
u.x[left]   = u.x[];
u.y[top]    = u.y[];
u.y[bottom] = u.y[];

double terrain (double x, double y)
{
  x -= -0.25; y -= -0.25;
  return 0.5*exp(-200.*(x*x + y*y));
}

void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = terrain (x, y);
}

event init (i = 0)
{
  zb.refine = refine_zb; // updates terrain
  foreach() {
    h[] = exp(-200.*(x*x + y*y));
    zb[] = terrain (x, y);
    h[] = max(h[] - zb[], 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);
  //  assert (s.min >= 0.);
}

event outputfile (t <= 1.2; t += 1.2/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l}, stdout, N);
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
