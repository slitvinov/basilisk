// Time-reversed advection in a vortex

#include "advection.h"
#include "vof.h"

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;
int MAXLEVEL;

void parameters() {
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
}

#define circle(x,y) (sq(0.2) - (sq(x + 0.2) + sq(y + .236338)))
#define END 15.

event init (i = 0) {
  scalar phi[];
  foreach_vertex()
    phi[] = circle(x,y);
  staggered vector s[];
  fractions (phi, f, s);
}

event velocity (i++) {
#if QUADTREE
  adapt_wavelet ({f}, (double[]){1e-2}, MAXLEVEL, list = {f});
#endif

  scalar psi[];
  foreach_vertex()
    psi[] = - 1.5*sin(2.*pi*t/END)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;
  trash ({u});
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[])/Delta;
  boundary ((scalar *){u});
}

event logfile (t = {0,END}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event field (t = END) {
  scalar phi[], e[];
  foreach_vertex()
    phi[] = circle(x,y);
  staggered vector s[];
  fractions (phi, e, s);
  foreach()
    e[] -= f[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

event shape (t += END/4.) {
  if (N == 128)
    output_facets (f, stdout);
}

#if 0
event movie (i += 10)
{
  static FILE * fp = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp, 1 << MAXLEVEL);
}
#endif

int main() {
  for (MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++) {
    N = 1 << MAXLEVEL;
    run();
  }
}
