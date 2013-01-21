/* definition of halo cells after coarsening */

#include <assert.h>

struct _Data {
  double h, e, w;
};

#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"
#include "adapt.c"

int main (int argc, char ** argv)
{
  void * grid = init_grid (32);

  double R0 = 0.1;
  foreach (grid) { data(0,0).h = exp(-(x*x + y*y)/(R0*R0)); } end_foreach();
  symmetry (grid, var(h));
  
  /* initial coarsening */
  restriction (grid, var(h));
  wavelet (grid, var(h), var(w));
  coarsen_wavelet (grid, var(w), 1e-2);
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
  } end_foreach_cell();

  free_grid (grid);
}
