#include "utils.h"

static int refine (Point point, void * data)
{
  return x*x + y*y < 0.25*0.25;
}

int main()
{
  init_grid(4);
  refine_function (refine, NULL, all);
  foreach_face(x)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face(y)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face()
    fprintf (stderr, "%g %g\n", x, y);
  output_cells (stdout);
  free_grid();
}
