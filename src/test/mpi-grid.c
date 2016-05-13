#include "grid/multigrid.h"

int main()
{
  foreach_dimension()
    periodic (right);

  L0 = 2.*pi;
  
  init_grid (4);
  output_cells (stdout);

  scalar s[];
  foreach()
    s[] = sin(x)*cos(y);
  boundary ({s});

  foreach()
    foreach_neighbor() {
      // fprintf (stderr, "%g %g %g %g\n", x, y, s[], s[] - sin(x)*cos(y));
      assert (fabs (s[] - sin(x)*cos(y)) < 1e-14);
    }
}
