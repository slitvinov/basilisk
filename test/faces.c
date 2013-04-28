#include "grid/quadtree.h"
#include "utils.h"

int main()
{
  init_grid(8);
  foreach_face(x)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face(y)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face()
    fprintf (stderr, "%g %g\n", x, y);
  output_cells (stdout);
  free_grid();
}
