/* speed of definition of halos */

#include <assert.h>
#include <time.h>

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
  
  clock_t start, end0, end;
  start = clock ();
  int i;
  for (i = 0; i < 10; i++) {
    /* coarsening */
    restriction (m, n, var(h));
    wavelet (m, n, var(h), var(w));
    coarsen_wavelet (m, n, var(w), 1e-5);
    flag_halo_cells (m, n);
    if (i == 0)
      end0 = clock();
  }
  end = clock ();
  double cpu0 = ((double) (end0 - start))/CLOCKS_PER_SEC;
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRID " %d iterations\n", i);
  fprintf (stderr, "initial coarsening: %6g CPU, %.3g points.steps/s\n",
	   cpu0, (n*n/cpu0));
  int leaves = 0, maxlevel = 0;
  foreach (m, n) { leaves++; if (level > maxlevel) maxlevel = level; } end_foreach();
  fprintf (stderr, "after coarsening: %d leaves, maximum level %d\n", leaves, maxlevel);
  fprintf (stderr, "iterations:         %6g CPU, %.3g leaves.steps/s\n",
	   cpu - cpu0, (leaves*(i - 1)/(cpu - cpu0)));

  FILE * fp = fopen("/tmp/halo", "w");
  foreach_cell (m, n) {
    if (!(cell.flags & halo))
      continue;
    else if (cell.flags & inactive)
      fprintf (fp, "%g %g %d %d\n", x, y, level, cell.neighbors);
  } end_foreach_cell();
  fclose(fp);

  free_grid (m);
}
