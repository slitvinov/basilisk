#include "grid/multigrid.h"
#include "utils.h"

static int refine_func (Point point, void * data)
{
  return x*x + y*y < 0.25*0.25;
}

int main()
{
  X0 = Y0 = -0.5;
  init_grid(4);
  //  refine_function (refine_func, NULL, NULL);
  output_cells (stdout);
  foreach()
    fprintf (stderr, "a %g %g\n", x - delta/2., y - delta/2.);
  foreach_boundary_face_ghost(right)
    fprintf (stderr, "b %g %g\n", x, y - delta/2.);
  foreach_boundary_face_ghost(top)
    fprintf (stderr, "c %g %g\n", x - delta/2., y);
  free_grid();
}
