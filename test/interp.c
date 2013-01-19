/* interpolation on halos  */

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
  int n = 1024;
  void * m = init_grid (n);

  double R0 = 0.1;
  foreach (m, n) { h(0,0) = exp(-(x*x + y*y)/(R0*R0)); } end_foreach();
  symmetry (m, n, var(h));
  
  /* initial coarsening (see halo.c) */
  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  coarsen_wavelet (m, n, var(w), 1e-4);
  flag_halo_cells (m, n);
  update_halos (m, n, var(h), var(h));

  FILE 
    * fp = fopen("/tmp/interp", "w");
  foreach_cell (m, n) {
    if (!(cell.flags & halo))
      continue;
    else if (cell.flags & inactive)
      fprintf (fp, "%g %g %d %d %g\n", x, y, level, cell.neighbors, h(0,0));
  } end_foreach_cell();

  fclose (fp);
  free_grid (m);
}
