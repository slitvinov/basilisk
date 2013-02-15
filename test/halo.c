/* definition of halo cells after coarsening */

#include <assert.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

new var h, w;

int main (int argc, char ** argv)
{
  void * grid = init_grid (32);

  double R0 = 0.1;
  foreach (grid)
    h(0,0) = exp(-(x*x + y*y)/(R0*R0));
  symmetry (grid, h);
  
  /* initial coarsening */
  restriction (grid, h);
  wavelet (grid, h, w);
  coarsen_wavelet (grid, w, 1e-2);
  flag_halo_cells (grid);

  foreach_cell (grid) {
    fprintf (stderr, "%g %g %d %d traversed\n", x, y, level, cell.neighbors);
    printf ("%g %g %d %d 1\n", x, y, level, cell.neighbors);
    if (!(cell.flags & halo))
      continue;
    else if (!(cell.flags & active)) {
      fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
      printf ("%g %g %d %d 2\n", x, y, level, cell.neighbors);
    }
    else {
      fprintf (stderr, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
      printf ("%g %g %d %d 3\n", x, y, level, cell.neighbors);
    }
  }

  free_grid (grid);
}
