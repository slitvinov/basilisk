/* definition of halo cells after refinement */

#include <assert.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

var h = new var, w = new var;

void refineiter (void * grid)
{
  for (int n = 0; n < 2; n++) {
    fprintf (stderr, "\nwavelet refinement\n");
    foreach(grid)
      h[] = exp(-(x*x + y*y)/(0.01));
    symmetry (grid, h);
    update_halo (grid, -1, h, h);

    restriction (grid, h, h);
    wavelet (grid, h, w);

    int nf = refine_wavelet (grid, h, h, w, 1e-2);
    flag_halo_cells (grid);

    fprintf (stderr, "refined %d cells\n", nf);
  }
  update_halo (grid, -1, h, h);
}

int main (int argc, char ** argv)
{
  void * grid = init_grid (16);

  refineiter (grid);

  foreach_halo(grid)
    printf ("%g %g %d %d %g %g halo1\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf(grid)
    printf ("%g %g %d %d %g %g leaf1\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  restriction (grid, h, h);
  wavelet (grid, h, w);
  fprintf (stderr, "\ncoarsened %d cells back\n", coarsen_wavelet (grid, w, 1e-2));
  flag_halo_cells (grid);

  foreach_halo(grid)
    printf ("%g %g %d %d %g %g halo2\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf(grid)
    printf ("%g %g %d %d %g %g leaf2\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  refineiter (grid);

  foreach_cell (grid) {
    if (!(cell.flags & halo))
      continue;
    else if (!(cell.flags & active))
      fprintf (stderr, "%g %g %d %d halo4\n", x, y, level, cell.neighbors);
    else
      fprintf (stderr, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
  }

  foreach_halo(grid)
    printf ("%g %g %d %d %g %g halo3\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf(grid)
    printf ("%g %g %d %d %g %g leaf3\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  free_grid (grid);
}
