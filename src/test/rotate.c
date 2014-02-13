// #include "grid/cartesian.h"
#include "advection.h"
#include "vof.h"

scalar c[];
scalar * interfaces = {c}, * tracers = NULL;
int MAXLEVEL;

int main()
{
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  for (MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++) {
    N = 1 << MAXLEVEL;
    run ();
  }
}

#define circle(x,y) (sq(0.1) - (sq(x-0.25) + sq(y)))

event init (i = 0)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = circle(x,y);
  fractions (phi, c);
}

#define end 0.785398

event velocity (i++) {
#if QUADTREE
  double cmax = 1e-2;
  adapt_wavelet ({c}, &cmax, MAXLEVEL, list = {c});
#endif

  trash ({u});
  foreach_face(x) u.x[] = -8.*y;
  foreach_face(y) u.y[] =  8.*x;
  boundary ((scalar *){u});
}

event logfile (t = {0,end}) {
  stats s = statsf (c);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event interface (t += end/10.) {
  static FILE * fp = fopen ("interface", "w");
  if (N == 64) {
    output_facets (c, fp);
    if (t == end)
      output_cells (fp);
  }
}

event field (t = end) {
  scalar e[], phi[];
  foreach_vertex()
    phi[] = circle(x,y);
  fractions (phi, e);
  foreach()
    e[] -= c[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  if (N == 64)
    output_field ({e}, stdout, N, linear = false);
}
