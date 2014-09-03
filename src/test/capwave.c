#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

scalar c[], * interfaces = {c}, * tracers = NULL;

uf.x[left]   = dirichlet(0);
uf.x[right]  = dirichlet(0);
uf.y[top]    = dirichlet(0);
uf.y[bottom] = dirichlet(0);

double se = 0, ne = 0;

int main() {
  L0 = 2.;
  Y0 = -L0/2.;
  c.sigma = 1.;
  TOLERANCE = 1e-6;
  const face vector muc[] = {0.0182571749236, 0.0182571749236};
  mu = muc;
  for (N = 16; N <= 256; N *= 2) {
    se = ne = 0;
    run();
  }
}

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - 0.01*cos (2.*pi*x));
  fractions (phi, c);
}

vector h[];

event amplitude (t += 3.04290519077e-3; t <= 2.2426211256) {
  heights (c, h);
  double max = - HUGE;;
  foreach() 
    if (c[] > 0 && c[] < 1) {
      double yi = y + h.y[];
      if (yi > max)
	max = yi;
    }
  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  static FILE * fp1 = fopen ("../prosperetti", "r");
  double t1, max1;
  fscanf (fp1, "%lf %lf", &t1, &max1);
  se += sq(max - max1); ne++;
  fprintf (fp, "%g %g\n", t*11.1366559937, max);
  fflush (fp);
}

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0, sqrt(se/ne)/0.01);

#if QUADTREE
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D -s oscillation.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set output 'amplitude.png'
set xlabel 'tau'
set ylabel 'Relative amplitude'
plot '../prosperetti' w l t "Prosperetti", 'wave-128' every 10 w p t "Basilisk"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set output 'convergence.png'
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
set grid
plot [5:200][1e-4:1]'clog' t "Basilisk" w lp, 2./x**2 t "Second order"
~~~
*/
