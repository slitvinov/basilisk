#include "grid/octree.h"

#if dimension == 2
double z = 0.;
#endif

int main()
{
  init_grid (8);
  origin (-0.5, -0.5, -0.5);
  mask (sq(x) + sq(y) + sq(z) < sq(0.4) ? none : right);
  output_cells (stdout);

  scalar a[];
  vector v[];
  face vector f[];

  a[right] = dirichlet(x + y + z);
  f.n[right] = 0.;
  
  foreach() {
    a[] = x + y + z;
    foreach_dimension()
      v.x[] = 1.;
  }

  boundary ({a});

  output_cells (stdout);
  foreach_boundary (depth()) {
    fprintf (stderr, "%g %g %g %g %g %g\n", x, y, z, a[], v.x[], v.y[]);
    assert (a[] == x + y + z);
  }
}
