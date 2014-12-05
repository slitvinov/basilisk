#include "fractions.h"
#include "heights.h"
#include "utils.h"

int main()
{
  init_grid (16);
  X0 = -0.5 - 0.125;
  Y0 = -0.5 - 0.125;
  scalar c[];
#if QUADTREE
  refine (level == 4 && fabs (x) < 0.375 && fabs (y) < 0.375, NULL);
  refine (level <= 5 && fabs (x) < 0.31 && fabs (y) < 0.31, NULL);
  c.refine = c.prolongation = fraction_refine;
#endif
  vertex scalar phi[];
  foreach_vertex()
    phi[] = - (0.2 - sqrt(sq(x+0.2) + sq(y+0.2)));
  fractions (phi, c);
  output_cells (stdout);
  output_facets (c, stdout);
  scalar kappa[];
  curvature (c, kappa);
  foreach()
    if (kappa[] != nodata)
      fprintf (stderr, "%g %g %g\n", x, y, kappa[]);
#if 0
  FILE * fp = popen ("gfsview2D -s", "w");
  output_gfs (fp);
#endif
}
