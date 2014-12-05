/* interpolation */

#include "utils.h"

scalar v[];

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

  refine (level == 3 && sq(x) + sq(y) > sq(0.49), NULL);
  refine (level <= 4 && sq(x) + sq(y) > sq(0.55), NULL);

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
