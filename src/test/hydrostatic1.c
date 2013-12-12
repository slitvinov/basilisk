// similar to hydrostatic.c but with a linear density profile
// (i.e. a quadratic pressure profile).

#include "navier-stokes/centered.h"

void parameters() {
  X0 = Y0 = -0.5;
  N = 8;
  DT = 1.;
  TOLERANCE = 1e-6;
  stokes = true;
}

p[top] = dirichlet(0);

int refine_func (Point point, void * data) {
  return (sq(x) + sq(y) < sq(0.25));
}

event init (i = 0) {
  refine_function (refine_func, NULL, {p,pf,u});
  const face vector g[] = {0,-1.};
  a = g;
  alpha = new face vector;
  foreach_face()
    alpha.x[] = 1./(0.51 - y);
}

event check (i = 1) {
  output_cells (stdout);
#if 1
  foreach()
    fprintf (stderr, "%g %g %g %.6f %.6f\n", x, y, p[]/dt, u.x[], u.y[]);
#else
  foreach_face(y)
    fprintf (stderr, "%g %g %g\n", x, y, uf.y[]);
#endif
}

int main() { run(); }
