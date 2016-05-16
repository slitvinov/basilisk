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
    s[] = sin(x)*cos(y)*cos(z);
  boundary ({s});

  foreach()
    foreach_neighbor(1) {
#if 0   
      fprintf (stderr, "%g %g %g %g %g\n", x, y, z, s[],
	       s[] - sin(x)*cos(y)*cos(z));
#else
      assert (fabs (s[] - sin(x)*cos(y)*cos(z)) < 1e-14);
#endif
    }
}
