#include "saint-venant.h"

#define LEVEL 7

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (i = 0)
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 2.5; t += 2.5/8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l});

#if !_MPI
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
#endif
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
