#include "terrain.h"
#include "utils.h"

int main ()
{
  FILE * fp = popen ("../kdt/xyz2kdt terrain", "w");
  for (double x = 0.; x <= 1.1; x += 0.005)
    for (double y = 0.; y <= 1.1; y += 0.005)
      fprintf (fp, "%g %g %g\n", x, y, sin(3.*pi*x)*cos(2.*pi*y));
  fclose (fp);

  for (int l = 4; l <= 7; l++) {
    init_grid (1 << l);
    scalar zb[];
    terrain (zb, "terrain");
    scalar e[];
    foreach()
      e[] = zb[] - sin(3.*pi*x)*cos(2.*pi*y);
    if (l == 7)
      output_field ({zb, e}, 128, stdout, false);
    norm n = normf (e);
    stats s = statsf (_terrain[zb].n);
    fprintf (stderr, "%d %g %g %g %g %g\n", l, n.avg, n.rms, n.max, 
	     s.min, s.max);
    free_grid();
  }
}