#include <assert.h>

struct _Data {
  double h, e, w;
};
#define h(k,l) data(k,l).h

#include "utils.h"
#include "quadtree.c"
#include "utils.c"

void restriction (Data * m, int n, var a)
{
  foreach_fine_to_coarse (m, n) {
    coarse(a,0,0) = (fine(a,0,0) + fine(a,1,0) + fine(a,0,1) + fine(a,1,1))/4.;
  } end_foreach_fine_to_coarse();
}

void wavelet (Data * m, int n, var a, var w)
{
  foreach_fine_to_coarse (m, n) {
    /* difference between fine value and bilinearly-interpolated coarse value */
    fine(w,0,0) = fine(a,0,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,-1)) + coarse(a,-1,-1))/16.;
    fine(w,0,1) = fine(a,0,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,+1)) + coarse(a,-1,+1))/16.;
    fine(w,1,0) = fine(a,1,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,0,-1) + coarse(a,+1,0)) + coarse(a,+1,-1))/16.;
    fine(w,1,1) = fine(a,1,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,+1,0) + coarse(a,0,+1)) + coarse(a,+1,+1))/16.;
  } end_foreach_fine_to_coarse ();
}

void coarsen_wavelet (Data * m, int n, var w, double max)
{
  foreach_fine_to_coarse (m, n) {
    double error = 0.;
    foreach_fine(k,l) {
      double e = fabs(fine(w,k,l));
      if (e > error)
	error = e;
    }
    if (error < max) {
      /* coarsen */
      cell.flags |= leaf;
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++) {
	  finecell(k,l).flags &= !leaf;
	  finecell(k,l).flags |= inactive;
	  /* update fine neighborhood */
	  for (int o = -GHOSTS; o <= GHOSTS; o++)
	    for (int p = -GHOSTS; p <= GHOSTS; p++)
	      finecell(k+o,l+p).neighbors--;
	}
    }
    /* propagate the error to coarser levels */
    coarse(w,0,0) = fabs(coarse(w,0,0)) + error;
  } end_foreach_fine_to_coarse ();
}

void flag_halo_cells (Data * m, int n)
{
  /* from the bottom up */
  foreach_cell_post (m, n, !(cell.flags & inactive) || cell.neighbors > 0) {
    if (cell.flags & inactive) {
      /* inactive and neighbors > 0 => this is a halo cell */
      cell.flags |= halo;
      /* propagate to parent */
      parentcell.flags |= halo;
    }
    else if (cell.flags & halo)
      /* propagate to parent */
      parentcell.flags |= halo;
  } end_foreach_cell_post();
}

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
