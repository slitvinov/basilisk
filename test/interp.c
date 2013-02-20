/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

var h = new var, w = new var;

int main (int argc, char ** argv)
{
  int n = 2048;
  void * grid = init_grid (n);

  double R0 = 0.1;
  foreach (grid)
    h[] = exp(-(x*x + y*y)/(R0*R0));
  symmetry (grid, h);
  
  /* initial coarsening (see halo.c) */
  restriction (grid, h);
  wavelet (grid, h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (grid, w, tolerance);
  flag_halo_cells (grid);

  restriction (grid, h);
  update_halo (grid, -1, h, h);

  double max = 0.;
  foreach_halo(grid) {
    double e = exp(-(x*x+y*y)/(R0*R0)) - h(0,0);
    printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h(0,0), e);
    if (fabs(e) > max)
      max = fabs(e);
  }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  free_grid (grid);

  return (max > tolerance);
}
