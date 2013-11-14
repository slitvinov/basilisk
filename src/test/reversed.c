// Time-reversed advection in a vortex

#include "grid/cartesian.h"
#include "advection.h"
#include "vof.h"

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

void parameters() {
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
}

#define circle(x,y) ((sq(x + 0.2) + sq(y + .236338)) < sq(0.2))

event init (i = 0) {
  foreach()
    f[] = circle(x,y);
  boundary ({f});
}

event velocity (i++) {
  trash ({u});
  foreach_face(x)
    u.x[] = 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*cos((y + 0.5)*pi);
  foreach_face(y)
    u.y[] = - 1.5*sin(2.*pi*t/5.)*cos((x + 0.5)*pi)*sin((y + 0.5)*pi);
  boundary ((scalar *){u});
}

event logfile (t = {0,5}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event field (t = 5) {
  scalar e[];
  foreach()
    e[] = f[] - circle(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

event shape (t += 0.5) {
  if (N == 256)
    output_facets (f, stdout);
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run();
}
