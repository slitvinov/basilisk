/* definition of halo cells after coarsening */

#include <assert.h>
#include "grid/quadtree.h"
#include "utils.h"

scalar h = new scalar, w = new scalar;

int main (int argc, char ** argv)
{
  init_grid (32);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary (h);
  
  /* initial coarsening */
  restriction (scalars(h));
  wavelet (h, w);
  coarsen_wavelet (w, 1e-2, 0);

  foreach_halo () {
    fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
    printf ("%g %g %d %d 2\n", x, y, level, cell.neighbors);
  }

  free_grid ();
}
