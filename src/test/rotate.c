#include "grid/cartesian.h"
#include "advection.h"
#include "vof.h"

scalar c[];
scalar * interfaces = {c}, * tracers = NULL;

void parameters()
{
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
}

#define circle(x,y) (sq(0.1) - (sq(x-0.25) + sq(y)))

event init (i = 0)
{
  scalar phi[];
  foreach_vertex()
    phi[] = circle(x,y);
  staggered vector s[];
  fractions (phi, c, s);
}

#define end 0.785398

event velocity (i++) {
  trash ({u});
  foreach_face(x) u.x[] = -8.*y;
  foreach_face(y) u.y[] =  8.*x;
  boundary ((scalar *){u});
}

event logfile (t = {0,end}) {
  stats s = statsf (c);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

#if 0
event interface (t += end/10.) {
  if (N == 256)
    output_facets (c, stdout);
}
#endif

event field (t = end) {
  scalar e[], phi[];
  foreach_vertex()
    phi[] = circle(x,y);
  staggered vector s[];
  fractions (phi, e, s);
  foreach()
    e[] -= c[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
#if 1
  if (N == 256)
    output_field ({e}, stdout, N, linear = false);
#endif
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run ();
}
