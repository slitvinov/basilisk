/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

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
  restriction (h, h);
  wavelet (h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (w, tolerance);
  flag_halo_cells (grid);

  restriction (h, h);
  update_halo (-1, h, h);

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
