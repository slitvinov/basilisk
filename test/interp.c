/* interpolation on halos  */

#include <assert.h>
#include <time.h>

struct _Data {
  double h, w;
};

#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

#define h(k,l) data(k,l).h

int main (int argc, char ** argv)
{
  int n = 2048;
  void * grid = init_grid (n);

  double R0 = 0.1;
  foreach (grid) { h(0,0) = exp(-(x*x + y*y)/(R0*R0)); }
  symmetry (grid, var(h));
  
  /* initial coarsening (see halo.c) */
  restriction (grid, var(h));
  wavelet (grid, var(h), var(w));
  double tolerance = 1e-4;
  coarsen_wavelet (grid, var(w), tolerance);
  flag_halo_cells (grid);

  update_halos (grid, var(h), var(h));

  FILE 
    * fp = fopen("/tmp/interp", "w");

  double max = 0.;
  foreach_halo(grid) {
    double e = exp(-(x*x+y*y)/(R0*R0)) - h(0,0);
    fprintf (fp, "%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h(0,0), e);
    if (fabs(e) > max)
      max = fabs(e);
  }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  fclose (fp);
  free_grid (grid);

  return (max > tolerance);
}
