/* definition of halo cells after coarsening */

#include "grid/quadtree.h"
#include "utils.h"

scalar h[];

int main (int argc, char ** argv)
{
  init_grid (32);

  X0 = Y0 = -0.5;
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
  
  /* initial coarsening */
  
  adapt_wavelet ({h}, (double []){1e-2}, 5);

  foreach_halo_coarse_to_fine(-1)
    fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
  foreach_halo_fine_to_coarse()
    fprintf (stderr, "%g %g %d %d res\n", x, y, level, cell.neighbors);
  for (int d = 0; d < nboundary; d++)
    foreach_boundary (d, true)
      fprintf (stderr, "%g %g %d %d boundary\n", x, y, level, cell.neighbors);
  output_cells (stdout);

  free_grid ();
}
