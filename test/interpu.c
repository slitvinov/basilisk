/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"
#include "wavelet.h"
#include "adapt.h"

scalar h = new scalar, u = new scalar, v = new scalar, w = new scalar;

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  symmetry (h);
  
  /* initial coarsening (see halo.c) */
  restriction (h, h);
  wavelet (h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (w, tolerance);
  flag_halo_cells (grid);

  foreach() {
    u[] = exp(-(xu*xu + yu*yu)/(R0*R0));
    v[] = exp(-(xv*xv + yv*yv)/(R0*R0));
  }

  restriction_u_v (u, v);
  update_halo_u_v (-1, u, v);

  double max = 0.;
  foreach() {
    double e = exp(-(x*x+y*y)/(R0*R0)) - (u[] + u[1,0])/2.;
    if (fabs(e) > max)
      max = fabs(e);
    printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, (u[] + u[1,0])/2., e);
  }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  free_grid ();

  return (max > 2.*tolerance);
}
