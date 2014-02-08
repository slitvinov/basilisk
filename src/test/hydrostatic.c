#include "navier-stokes/centered.h"

int main() {
  X0 = Y0 = -0.5;
  N = 8;
  DT = 1.;
  TOLERANCE = 1e-6;
  stokes = true;
  run();
}

int refine_func (Point point, void * data) {
  return (sq(x) + sq(y) < sq(0.25));
}

event init (i = 0) {
  refine_function (refine_func, NULL, {p,pf,u});
  const face vector g[] = {0,-1.};
  a = g;
}

event check (i = 1) {
  output_cells (stdout);
#if 1
  foreach()
    fprintf (stderr, "%g %g %g %.6f %.6f\n", x, y, p[],
	     fabs(u.x[]), fabs(u.y[]));
#else
  foreach_face(y)
    fprintf (stderr, "%g %g %g\n", x, y, uf.y[]);
#endif
  foreach()
    assert (fabs (u.y[]) < 1e-6);
}
