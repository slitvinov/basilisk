#include "grid/cartesian.h"
#include "fractions.h"

int main()
{
  X0 = Y0 = -0.5;
  init_grid (32);

  scalar c[], phi[];
  staggered vector s[];
  foreach_vertex()
    phi[] = sq(0.25) - sq(x) - sq(y);
  fractions (phi, c, s);

  output_fractions (c, s, stdout);
  output_facets (c, stderr);

  free_grid ();
}
