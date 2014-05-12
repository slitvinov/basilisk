/* interpolation */

#include "utils.h"

scalar v[];

double radius;

int refine_func (Point point, void * data) {
  return sq(x) + sq(y) > sq(radius);
}

void boundary_halo (int d, int l, int depth) {
  foreach_boundary_cell (d, true) {
    if (level == l) {
      if (level == depth ||
	  is_leaf(cell) || (cell.flags & halo) || is_corner(cell)) {
	fprintf (stderr, "%g %g %s\n", x, y,
		 level == depth ?  "D" :
		 is_leaf(cell) ?   "L" : 
		 is_corner(cell) ? "C" : 
		                   "R");
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	foreach_child_direction(d)
	  fprintf (stderr, "%g %g P\n", x, y);
      continue;
    }
  }
}

int main (int argc, char ** argv)
{
  origin (-0.5, -0.5);
  init_grid (8);

  radius = 0.49;
  refine_function (refine_func, NULL, NULL);
  radius = 0.55;
  refine_function (refine_func, NULL, NULL);

  output_cells (stdout);

  FILE * fp = fopen ("restriction", "w");
  for (int l = 2; l < 5; l++)
    foreach_halo (restriction, l)
      fprintf (fp, "%d %g %g\n", l, x, y);
  fclose (fp);

  fp = fopen ("prolongation", "w");
  for (int l = 2; l < 5; l++)
    foreach_halo (prolongation, l)
      fprintf (fp, "%d %g %g\n", l, x, y);
  fclose (fp);

  boundary_halo (top, 3, 4);
  boundary_halo (left, 4, 5);
  boundary_halo (right, 5, 6);

  free_grid();
}
