// similar to fractions.tst but with refinement
#include "fractions.h"

static int band (Point point, void * data)
{
  return fabs (x) < 0.25;
}

int main()
{
  X0 = Y0 = -0.5;
  init_grid (16);

  refine_function (band, NULL, all);

  scalar c[], phi[];
  staggered vector s[];
  foreach_vertex()
    phi[] = sq(0.3) - sq(x) - sq(y);
  fractions (phi, c, s);

  output_fractions (c, s, stdout);
  output_facets (c, stderr);
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);

  vector n[];
  scalar alpha[];
  reconstruction (c, n, alpha);
  coord p[2];
  foreach_halo() {
    coord m = {n.x[],n.y[]};
    if (facets (c[], m, alpha[], p) == 2)
      fprintf (stderr, "halo %g %g\nhalo %g %g\nhalo\n", 
	       x + p[0].x*Delta, y + p[0].y*Delta, 
	       x + p[1].x*Delta, y + p[1].y*Delta);
  }

  free_grid ();
}
