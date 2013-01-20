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
  int n = 32;
  void * m = init_grid (n);

  double R0 = 0.1;
  foreach (m, n) { data(0,0).h = exp(-(x*x + y*y)/(R0*R0)); } end_foreach();
  symmetry (m, n, var(h));
  
  /* initial coarsening */
  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  coarsen_wavelet (m, n, var(w), 1e-2);
  flag_halo_cells (m, n);

  foreach_cell (m, n) {
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

  free_grid (m);
}
