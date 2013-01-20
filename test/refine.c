/* definition of halo cells after refinement */

#include <assert.h>

struct _Data {
  double h, w;
};

#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"
#include "adapt.c"

int refine (void * m, int n, var w, double max, var start, var end)
{
  int nf = 0;
  foreach_leaf (m, n) {
    if (fabs(val(w,0,0)) >= max) {
      alloc_children();
      cell.flags &= ~leaf;
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++) {
	  assert(!(child(k,l).flags & active));
	  child(k,l).flags |= (active | leaf);
	  /* bilinear interpolation from coarser level */
	  for (var a = start; a <= end; a += sizeof(double)) /* for each variable */
	    fine(a,k,l) = 
	      (9.*val(a,0,0) + 3.*(val(a,2*k-1,0) + val(a,0,2*l-1)) + val(a,2*k-1,2*l-1))/16.;
	  /* update neighborhood */
	  for (int o = -GHOSTS; o <= GHOSTS; o++)
	    for (int p = -GHOSTS; p <= GHOSTS; p++)
	      child(k+o,l+p).neighbors++;
	}
      nf++;
    }
  } end_foreach_leaf();
  return nf;
}

void refineiter (void * m, int n)
{
  for (int n = 0; n < 2; n++) {
    fprintf (stderr, "\nwavelet refinement\n");
    foreach(m,n) {
      data(0,0).h = exp(-(x*x + y*y)/(0.01));
    } end_foreach();
    symmetry (m, n, var(h));
    update_halos (m, n, var(h), var(h));

    restriction (m, n, var(h));
    wavelet (m, n, var(h), var(w));

    int nf = refine (m, n, var(w), 1e-2, var(h), var(h));
    flag_halo_cells (m, n);

    fprintf (stderr, "refined %d cells\n", nf);
  }
  update_halos (m, n, var(h), var(h));
}

int main (int argc, char ** argv)
{
  void * m = init_grid (1);
  int n = 1; /* not used anyway */

  fprintf (stderr, "refine uniformly to 4 levels\n");
  for (int l = 0; l < 4; l++) {
    int nf = refine (m, n, var(w), 0., var(h), var(h));
    fprintf (stderr, "refined %d cells\n", nf);
  }

  refineiter (m, n);

  foreach_halo(m,n) {
    printf ("%g %g %d %d %g %g halo1\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(m,n) {
    printf ("%g %g %d %d %g %g leaf1\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  fprintf (stderr, "\ncoarsened %d cells back\n", coarsen_wavelet (m, n, var(w), 1e-2));
  flag_halo_cells (m, n);

  foreach_halo(m,n) {
    printf ("%g %g %d %d %g %g halo2\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(m,n) {
    printf ("%g %g %d %d %g %g leaf2\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  refineiter (m, n);

  foreach_cell (m, n) {
    if (!(cell.flags & halo))
      continue;
    else if (!(cell.flags & active))
      fprintf (stderr, "%g %g %d %d halo4\n", x, y, level, cell.neighbors);
    else
      fprintf (stderr, "%g %g %d %d flagged\n", x, y, level, cell.neighbors);
  } end_foreach_cell();

  foreach_halo(m,n) {
    printf ("%g %g %d %d %g %g halo3\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_halo();
  foreach_leaf(m,n) {
    printf ("%g %g %d %d %g %g leaf3\n", x, y, level, cell.neighbors, data(0,0).h,
	    fabs(data(0,0).h - exp(-(x*x + y*y)/(0.01))));
  } end_foreach_leaf();

  free_grid (m);
}
