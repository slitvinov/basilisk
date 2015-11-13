#include "grid/multigrid3D.h"
#include "fractions.h"
#include "curvature.h"
#include "utils.h"

int main()
{
  origin (-0.5, -0.5, -0.5);
  init_grid (16);
  scalar f[];
  vertex scalar phi[];
  foreach_vertex()
    phi[] = x*x + y*y + z*z - sq(0.25001);
  fractions (phi, f);
  vector h[];
  heights (f, h);
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  fprintf (stderr, "kappa min: %g avg: %g stddev: %g max: %g\n",
	   s.min, s.sum/s.volume, s.stddev, s.max);
#if 0  
  FILE * fp = popen ("gfsview3D -s hf.gfv", "w");
  output_gfs (fp);
#endif
}
