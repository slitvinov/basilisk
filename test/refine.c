/* definition of halo cells after refinement */

#include <assert.h>

struct _Data {
  double h, w;
};

#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"
#include "adapt.c"

void refineiter (void * grid)
{
  for (int n = 0; n < 2; n++) {
    fprintf (stderr, "\nwavelet refinement\n");
    foreach(grid) {
      data(0,0).h = exp(-(x*x + y*y)/(0.01));
    } end_foreach();
    symmetry (grid, var(h));
    update_halos (grid, var(h), var(h));

    restriction (grid, var(h));
    wavelet (grid, var(h), var(w));

    int nf = refine_wavelet (grid, var(w), 1e-2, var(h), var(h));
    flag_halo_cells (grid);

    fprintf (stderr, "refined %d cells\n", nf);
  }
  update_halos (grid, var(h), var(h));
}

int main (int argc, char ** argv)
{
  void * grid = init_grid (16);

  refineiter (grid);

  foreach_halo(grid) {
    printf ("%g %g %d %d %g %g halo1\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(grid) {
    printf ("%g %g %d %d %g %g leaf1\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  restriction (grid, var(h));
  wavelet (grid, var(h), var(w));
  fprintf (stderr, "\ncoarsened %d cells back\n", coarsen_wavelet (grid, var(w), 1e-2));
  flag_halo_cells (grid);

  foreach_halo(grid) {
    printf ("%g %g %d %d %g %g halo2\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(grid) {
    printf ("%g %g %d %d %g %g leaf2\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  refineiter (grid);

  foreach_cell (grid) {
    if (!(cell.flags & halo))
      continue;
    else if (!(cell.flags & active))
      fprintf (stderr, "%g %g %d %d halo4\n", x, y, level, cell.neighbors);
    else
      fprintf (stderr, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
  } end_foreach_cell();

  foreach_halo(grid) {
    printf ("%g %g %d %d %g %g halo3\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(grid) {
    printf ("%g %g %d %d %g %g leaf3\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  free_grid (grid);
}
