/* definition of halo cells after coarsening */

#include <assert.h>

struct _Data {
  double h, e, w;
};
#define h(k,l) data(k,l).h

#include "utils.h"
#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"
#include "adapt.c"

int main (int argc, char ** argv)
{
  int n = 32;
  void * m = init_grid (n);

  double R0 = 0.1;
  foreach (m, n) { h(0,0) = exp(-(x*x + y*y)/(R0*R0)); } end_foreach();
  symmetry (m, n, var(h));
  
  /* initial coarsening */
  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  coarsen_wavelet (m, n, var(w), 1e-2);
  flag_halo_cells (m, n);

  FILE 
    * ftraversed = fopen("/tmp/traversed", "w"), 
    * fhalo = fopen("/tmp/halo", "w"),
    * flagged = fopen("/tmp/flagged", "w");

  foreach_cell (m, n) {
    fprintf (stderr, "%g %g %d %d traversed\n", x, y, level, cell.neighbors);
    fprintf (ftraversed, "%g %g %d %d traversed\n", x, y, level, cell.neighbors);
    if (!(cell.flags & halo))
      continue;
    else if (cell.flags & inactive) {
      fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
      fprintf (fhalo, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
    }
    else {
      fprintf (flagged, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
      fprintf (stderr, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
    }
  } end_foreach_cell();

  fclose (ftraversed);
  fclose (fhalo);
  fclose (flagged);

  free_grid (m);
}
