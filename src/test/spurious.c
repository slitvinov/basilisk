#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

scalar c[], * interfaces = {c}, * tracers = NULL;

#define DIAMETER 0.8
#define MU sqrt(DIAMETER/LAPLACE)
#define TMAX (sq(DIAMETER)/MU)

int LEVEL;
double LAPLACE;
double DC = 0.;

FILE * fp = NULL;

int main() {
  TOLERANCE = 1e-6;
  stokes = true;
  c.sigma = 1;
  LEVEL = 5;
  N = 1 << LEVEL;
  for (LAPLACE = 120; LAPLACE <= 12000; LAPLACE *= 10)
    run();

  LAPLACE = 12000; DC = 1e-10;
  for (LEVEL = 3; LEVEL <= 7; LEVEL++) 
    if (LEVEL != 5) {
      N = 1 << LEVEL;
      run();
    }
}

scalar cref[], cn[];

event init (i = 0) {
  const face vector muc[] = {MU,MU};
  mu = muc;

  char name[80];
  sprintf (name, "La-%g-%d", LAPLACE, LEVEL);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");
  
  scalar phi[];
  foreach_vertex()
    phi[] = sq(DIAMETER/2) - sq(x) - sq(y);
  fractions (phi, c);
  foreach()
    cref[] = c[];  
  boundary ({cref});

  foreach()
    cn[] = c[];
}

event logfile (i++; t <= TMAX)
{  
  double dc = change (c, cn);
  if (i > 0 && dc < DC)
    return 1; /* stop */

  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g\n", MU*t/sq(DIAMETER), normf(un).max*sqrt(DIAMETER));
}

event error (t = end) {
  double vol = statsf(c).sum;
  double radius = sqrt(4.*vol/pi);
  double ekmax = 0.;
  scalar un[], ec[];
  scalar kappa = c.kappa;
  foreach() {
    un[] = norm(u);
    ec[] = c[] - cref[];
    if (kappa[] != nodata) {
      double ek = fabs (kappa[] - 1./radius);
      if (ek > ekmax)
	ekmax = ek;
    }
  }
  norm ne = normf (ec);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", 
	   LEVEL, LAPLACE, 
	   normf(un).max*sqrt(DIAMETER), 
	   ne.avg, ne.rms, ne.max,
	   ekmax);
}

#if QUADTREE
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s spurious.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
~~~gnuplot Evolution of the amplitude of the capillary currents $\max(|{\bf u}|)(D/\sigma)^{1/2}$ as a function of non-dimensional time $\tau=t\mu/D^2$ for the range of Laplace numbers indicated in the legend.
set xlabel 't{/Symbol m}/D^2'
set ylabel 'U(D/{/Symbol s})^{1/2}'
set logscale y
plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", \
  'La-12000-5' w l t "La=12000"
~~~

~~~gnuplot Convergence of the error on the equilibrium shape of the droplet with resolution. The diameter is given in number of grid points.
set xlabel 'D'
set ylabel 'Shape error'
set logscale x
set xtics 2
plot [5:120]'< sort -n -k1,2 clog' u (0.8*2**$1):5 w lp t "RMS" ps 2, \
            '< sort -n -k1,2 clog' u (0.8*2**$1):6 w lp t "Max" ps 2, \
             0.2/(x*x) t "Second order"
~~~

~~~gnuplot Convergence of the relative error on the equilibrium curvature value with resolution. The diameter is given in number of grid points.
set ylabel 'Relative curvature error'
plot [5:120]'< sort -n -k1,2 clog' u (0.8*2**$1):($7/2.5) w lp t "Max" ps 2, \
             0.6/(x*x) t "Second order"
~~~
*/
