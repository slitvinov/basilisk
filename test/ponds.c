// Bumps uniformly covered with water create ponds and streams (i.e. a
// lot of wetting/drying going on)

// use KDT database rather than analytical function
#define KDT 1

#include "saint-venant.h"
#if KDT
# include "terrain.h"
#endif

#define LEVEL 7

void parameters()
{
  N = 1 << LEVEL;
  G = 9.81;
  L0 = 1000.;
}

#define zb(x,y) ((cos(pi*x/L0)*cos(pi*y/L0) + \
		  cos(3.*pi*x/L0)*cos(3.*pi*y/L0)) - 2.*x/1000.)

#if !KDT
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = zb(x,y);
}
#endif

void init()
{
#if !KDT
  zb.refine = refine_zb; // updates terrain
  foreach() {
    zb[] = zb(x,y);
    h[] = 0.1;
  }
#else
  FILE * fp = popen ("../kdt/xyz2kdt ponds", "w");
  for (double x = 0.; x <= 1000; x += 1.)
    for (double y = 0.; y <= 1000.; y += 1.)
      fprintf (fp, "%g %g %g\n", x, y, zb(x,y));
  fclose (fp);
  terrain (zb, "ponds", NULL);
  foreach()
    h[] = 0.1;
#endif
}

event friction (i++) {
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
  boundary ((scalar *){u});
}

event logfile (i += 10) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g %g\n", 
	   t, i, s.min, s.max, s.sum, dt, n.rms, n.max, dt);
  assert (s.min > 0.);
}

event outputfile (t <= 1200.; t += 1200./8) {
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h, zb, u}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({h, zb, l}, stdout, N);

  nf++;
}

// int event (t += 100)
//  output_matrix (h, stdout, N, true);

event adapt (i++) {
  scalar w[];
  wavelet (h, w);

  double cmax = 1e-2;
  int nf = refine_wavelet (w, cmax, LEVEL, all);
  int nc = coarsen_wavelet (w, cmax/4., 4, all);
  if (nf || nc)
    boundary (all);

  //  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}

int main() { run(); }
