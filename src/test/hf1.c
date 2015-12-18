#include "fractions.h"
#include "curvature.h"
#include "utils.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (16);
  scalar c[];
  vertex scalar phi[];
  foreach_vertex()
    phi[] = min(sq(x + 0.1) + sq(y + 0.22) - sq(0.18),
		sq(x - 0.1) + sq(y - 0.22) - sq(0.18));
  fractions (phi, c);
  output_cells (stdout);
  output_facets (c, stdout);
  vector h[];
  heights (c, h);
  scalar kappa[];
  curvature (c, kappa);
  foreach() {
    if (h.x[] != nodata)
      fprintf (stderr, "%g %g %g %g hx\n", x, y, height(h.x[]), Delta);
    if (h.y[] != nodata)
      fprintf (stderr, "%g %g %g %g hy\n", x, y, height(h.y[]), Delta);
    if (kappa[] != nodata)
      fprintf (stderr, "%g %g %g kappa\n", x, y, kappa[]);
  }
}