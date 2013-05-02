/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"

scalar h = new scalar, w = new scalar;

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary (h);
  
  /* initial coarsening (see halo.c) */
  restriction (h);
  boundary_all (scalars (h));
  wavelet (h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (w, tolerance, 0, none);

  halo_restriction (scalars (h));
  halo_prolongation (-1, scalars (h));

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
