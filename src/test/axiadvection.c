#include "navier-stokes/centered.h"
#include "axi.h"
#include "vof.h"
#include "tracer.h"

scalar f[], f1[];
scalar * interfaces = {f}, * tracers = {f1};

int main()
{
  X0 = -0.5;
  N = 64;
  TOLERANCE = 1e-6;
  f1.gradient = minmod2;
  run();
}

u.x[left] = dirichlet(1);
u.y[left] = dirichlet(0);
p[left]   = neumann(0);

u.y[top] = neumann(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);

#define ellipse(xc, yc, a, b) (sq((x - xc)/(a)) + sq((y - yc)/(b)) - 1.)

event init (i = 0) {
  foreach()
    u.x[] = 1.;
  scalar phi[];
  foreach_vertex()
    phi[] = - ellipse (0, 0.3, 0.1, 0.1);
  fractions (phi, f);
  fractions (phi, f1);
}

event logfile (i++; t <= 0.8) {
  static double sfmin = HUGE, sfmax = -HUGE;
  static double sfmin1 = HUGE, sfmax1 = -HUGE;
  double s = statsf(f).sum, s1 = statsf(f1).sum;
  if (s < sfmin) sfmin = s;
  if (s > sfmax) sfmax = s;
  if (s1 < sfmin1) sfmin1 = s1;
  if (s1 > sfmax1) sfmax1 = s1;
  double e = 2.*(sfmax - sfmin)/(sfmax + sfmin);
  double e1 = 2.*(sfmax1 - sfmin1)/(sfmax1 + sfmin1);
  fprintf (stderr, "%g %.12f %.12f %g %g\n", t, s, s1, e, e1);
  assert (e < 2e-6);
  assert (e1 < 5e-5);
}

event output (t += 0.2; t <= 1.2)
  output_facets (f);

#if QUADTREE

#if 0
event gfsview (i++) {
  static FILE * fp = popen ("gfsview2D -s test.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

event adapt (i++) {
  double sb = statsf(f).sum;
  double sb1 = statsf(f1).sum;
  /* we need to adapt p to get an initial guess for the next
     iteration, but we can't apply boundary conditions (for the
     pressure) because alpha is not consistent and is used to compute
     the consistent Neumann conditions (see navier-stokes/centered.h). */
  adapt_wavelet ({f1}, (double[]){5e-3}, 
		 maxlevel = 6, minlevel = 4,
		 list = {cm,fm,p,u,pf,uf,g,f,f1}, 
		 listb = {cm,fm,u,pf,uf,g,f,f1});
  double sa = statsf(f).sum;
  double sa1 = statsf(f1).sum;
  // the mass of VOF tracers is not conserved exactly
  assert (fabs(sa - sb) < 1e-6);
  // the mass of diffusive tracers must be conserved to within round-off
  assert (fabs(sa1 - sb1) < 1e-12);
}
#endif

/**
~~~gnuplot Evolution of the VOF interface
set size ratio -1
plot [-0.5:0.5][0:1]'out' w l t ''
~~~
*/
