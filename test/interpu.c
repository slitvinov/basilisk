/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

var h = new var, u = new var, v = new var, w = new var;

int main (int argc, char ** argv)
{
  int n = 2048;
  void * grid = init_grid (n);

  double R0 = 0.1;
  foreach (grid)
    h[] = exp(-(x*x + y*y)/(R0*R0));
  symmetry (grid, h);
  
  /* initial coarsening (see halo.c) */
  restriction (grid, h, h);
  wavelet (grid, h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (grid, w, tolerance);
  flag_halo_cells (grid);

  foreach (grid) {
    u[] = exp(-(xu*xu + yu*yu)/(R0*R0));
    v[] = exp(-(xv*xv + yv*yv)/(R0*R0));
  }

  restriction_u_v (grid, u, v);
  update_halo_u_v (grid, -1, u, v);

  double max = 0.;
  foreach (grid) {
    double e = exp(-(x*x+y*y)/(R0*R0)) - (u[] + u[1,0])/2.;
    if (fabs(e) > max)
      max = fabs(e);
    printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, (u[] + u[1,0])/2., e);
  }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  free_grid (grid);

  return (max > 2.*tolerance);
}
