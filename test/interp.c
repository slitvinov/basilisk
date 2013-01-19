/* interpolation on halos  */

#include <assert.h>
#include <time.h>

struct _Data {
  double h, w;
};

#include "utils.h"
#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"
#include "adapt.c"

#define h(k,l) data(k,l).h

int main (int argc, char ** argv)
{
  int n = 2048;
  void * m = init_grid (n);

  double R0 = 0.1;
  foreach (m, n) { h(0,0) = exp(-(x*x + y*y)/(R0*R0)); } end_foreach();
  symmetry (m, n, var(h));
  
  /* initial coarsening (see halo.c) */
  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  double tolerance = 1e-4;
  coarsen_wavelet (m, n, var(w), tolerance);
  flag_halo_cells (m, n);

  update_halos (m, n, var(h), var(h));

  FILE 
    * fp = fopen("/tmp/interp", "w");

  double max = 0.;
  foreach_halo(m, n) {
    double e = exp(-(x*x+y*y)/(R0*R0)) - h(0,0);
    fprintf (fp, "%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h(0,0), e);
    if (fabs(e) > max)
      max = fabs(e);
  } end_foreach_halo();

  fprintf (stderr, "maximum error on halos: %g\n", max);

  fclose (fp);
  free_grid (m);

  return (max > tolerance);
}
