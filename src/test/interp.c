/* interpolation on halos  */

#include "grid/quadtree.h"
#include "utils.h"

scalar h[];

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
  
  /* initial coarsening (see halo.c) */
  double tolerance = 1e-4;
  adapt_wavelet ({h}, &tolerance, 11);

  halo_restriction ({h}, NULL);
  halo_prolongation (-1, {h});

  double max = 0.;
  foreach_halo() {
    double e = exp(-(x*x+y*y)/(R0*R0)) - h[];
    printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h[], e);
    if (fabs(e) > max)
      max = fabs(e);
  }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  free_grid ();

  return (max > tolerance);
}
