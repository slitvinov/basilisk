/* speed of definition of halos */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"

scalar h = new scalar, w = new scalar;

int main (int argc, char ** argv)
{
  int n = 1024;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary (h);
  
  clock_t start, end0, end;
  start = end0 = clock ();
  int i;
  for (i = 0; i < 31; i++) {
    /* coarsening */
    restriction (scalars (h));
    wavelet (h, w);
    coarsen_wavelet (w, 1e-5, 0);
    if (i == 0)
      end0 = clock();
  }
  end = clock ();
  double cpu0 = ((double) (end0 - start))/CLOCKS_PER_SEC;
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "---- restriction + wavelet + coarsen_wavelet + flag_halo_cells ----\n");
  int leaves = 0, maxlevel = 0;
  foreach() { leaves++; if (level > maxlevel) maxlevel = level; }
  fprintf (stderr, "after coarsening: %d leaves, maximum level %d\n", leaves, maxlevel);
  fprintf (stderr, "initial coarsening:  %6g CPU, %.3g points.steps/s\n",
	   cpu0, n*n/cpu0);
  fprintf (stderr, "%4d iterations:     %6g CPU, %.3g leaves.steps/s\n",
	   i - 1, cpu - cpu0, leaves*(i - 1)/(cpu - cpu0));

  int nhalos = 0;
  foreach_halo() {
    printf ("%g %g %d %d\n", x, y, level, cell.neighbors);
    nhalos++;
  }

  start = clock ();
  for (i = 0; i < 200; i++)
    halo_interpolation (-1, scalars(h));
  end = clock ();
  cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "---- update_halos ----\n");
  fprintf (stderr, "%d halo points\n", nhalos);
  fprintf (stderr, "%4d iterations:     %6g CPU, %.3g halos.steps/s\n",
	   i, cpu, nhalos*i/cpu);

  free_grid();
}
