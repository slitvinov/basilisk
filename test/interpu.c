/* interpolation on halos  */

#include <assert.h>
#include <time.h>
#include "grid/quadtree.h"
#include "utils.h"

scalar h = new scalar, u = new scalar, v = new scalar, w = new scalar;

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary (h);
  
  /* initial coarsening (see halo.c) */
  wavelet (h, w);
  double tolerance = 1e-4;
  coarsen_wavelet (w, tolerance, 0, none);

  foreach_face(x) u[] = exp(-(x*x + y*y)/(R0*R0));
  foreach_face(y) v[] = exp(-(x*x + y*y)/(R0*R0));

  vector uv = {u,v};
  /* see boundary_quadtree() */
  halo_restriction_flux (vectors (uv));
  scalar * list = scalars (uv);
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_cell (b, true) {
      if (is_active (cell)) {
	if (cell.neighbors > 0)
	  for (scalar s in list)
	    s[ghost] = _boundary[b][s] (point, s);
      }
      else
	continue;
    }
  halo_prolongation_u_v (-1, u, v);

  double max = 0., maxv = 0;
  foreach_halo() {
    double xu = x - delta/2., yu = y;
    double xv = x, yv = y - delta/2.;
    double e = exp(-(xu*xu+yu*yu)/(R0*R0)) - u[];
    if (fabs(e) > max)
      max = fabs(e);
    printf ("%g %g %d %d %g %g\n", xu, yu, level, cell.neighbors, u[], e);
    e = fabs (exp(-(xv*xv+yv*yv)/(R0*R0)) - v[]);
    if (e > maxv)
      maxv = e;
  }

  fprintf (stderr, "maximum error on halos: %g %g\n", max, maxv);

  free_grid ();

  return (max > tolerance || maxv > tolerance || max != maxv);
}
