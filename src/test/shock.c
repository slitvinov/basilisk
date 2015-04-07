/* similar to:
   http://gerris.dalembert.upmc.fr/gerris/examples/examples/shock.html */

#include "saint-venant.h"

int LEVEL = 9;

bid cylinder;

int main()
{
  size (5.);
  G = 9.81;
  origin (-L0/2., -L0/2.);
  init_grid (1 << LEVEL);
  run();
}

#define H0 3.505271526
#define U0 6.29033769408481

h[left]   = H0;
eta[left] = H0;
u.n[left] = U0;

event init (i = 0)
{
  mask (sq(x + 0.5) + sq(y) < sq(0.5) ? cylinder : none);
  foreach() {
    h[] = (x <= -1 ? H0 : 1.);
    u.x[] = (x <= -1 ? U0 : 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event movie (t += 0.0025)
{
  static FILE * fp = popen ("ppm2mpeg > depth.mpg", "w");
  output_ppm (h, fp,
	      min = 0.1, max = 6, map = cool_warm, n = 400, linear = true);
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, map = cool_warm, min = 0, max = LEVEL, n = 400);
}

event image (t = end) {
  output_ppm (h, file = "depth.png",
	      min = 0.1, max = 6, map = cool_warm, n = 400, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "level.png", map = cool_warm,
	      min = 0, max = LEVEL, n = 400);
}
  
event outputfile (t <= 0.3; t += 0.3/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l});
}

#if 1
event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-2}, LEVEL, 5);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif
