/* definition of halo cells after refinement */

#include <assert.h>
#include "grid/quadtree.h"
#include "utils.h"

scalar h = new scalar, w = new scalar;

void refineiter ()
{
  for (int n = 0; n < 2; n++) {
    fprintf (stderr, "\nwavelet refinement\n");
    foreach()
      h[] = exp(-(x*x + y*y)/(0.01));
    boundary (h);

    restriction (h);
    wavelet (h, w);

    int nf = refine_wavelet (w, 1e-2, 12, scalars (h));

    fprintf (stderr, "refined %d cells\n", nf);
  }
  boundary (h);
}

int main (int argc, char ** argv)
{
  init_grid (16);

  refineiter();

  foreach_halo()
    printf ("%g %g %d %d %g %g halo1\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf()
    printf ("%g %g %d %d %g %g leaf1\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  restriction (h);
  wavelet (h, w);
  fprintf (stderr, "\ncoarsened %d cells back\n", 
	   coarsen_wavelet (w, 1e-2, 0, none));

  foreach_halo()
    printf ("%g %g %d %d %g %g halo2\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf()
    printf ("%g %g %d %d %g %g leaf2\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  refineiter();

  foreach_halo()
    fprintf (stderr, "%g %g %d %d halo4\n", x, y, level, cell.neighbors);

  foreach_halo()
    printf ("%g %g %d %d %g %g halo3\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));
  foreach_leaf()
    printf ("%g %g %d %d %g %g leaf3\n", x, y, level, cell.neighbors, h[],
	    fabs(h[] - exp(-(x*x + y*y)/(0.01))));

  free_grid ();
}
