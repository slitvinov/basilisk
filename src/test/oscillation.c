#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

scalar c[], * interfaces = {c}, * tracers = NULL;

#define D 0.2
#define rho(c) (1.*(c) + 1e-3*(1. - (c)))

FILE * fp = NULL;

int LEVEL;

int main() {
  L0 = 0.5;
  c.sigma = 1.;
  remove ("error");
  remove ("laplace");
  for (LEVEL = 5; LEVEL <= 7; LEVEL++) {
    N = 1 << LEVEL;
    char name[80];
    sprintf (name, "k-%d", LEVEL);
    fp = fopen (name, "w");
    run();
    fclose (fp);
  }
  system ("grep ^fit out >> log");
}

#if 0
#define cf c
#else
scalar cf[];
#endif

void density();

event init (i = 0) {
  alpha = new face vector;
  scalar phi[];
  foreach_vertex()
    phi[] = D/2.*(1. + 0.05*cos(2.*atan2(y,x))) - sqrt(sq(x) + sq(y));
  fractions (phi, c);
#ifndef cf
  foreach()
    cf[] = c[];
#endif
  density();
}

void density() {
#ifndef cf
  foreach()
    cf[] = (4.*c[] + 
	    2.*(c[0,1] + c[0,-1] + c[1,0] + c[-1,0]) +
	    c[-1,-1] + c[1,-1] + c[1,1] + c[-1,1])/16.;
  boundary ({cf});
#endif
  foreach_face() {
    double cm = (cf[] + cf[-1,0])/2.;
    alpha.x[] = 1./rho(cm);
  }
  boundary ((scalar *){alpha});
}

event properties (i++)
  density();

event logfile (i++; t <= 1) {
  double ke = 0.;
  foreach (reduction(+:ke))
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(cf[]);
  fprintf (fp, "%g %g %d\n", t, ke, mgp.i);
}

event fit (t = end) {
  FILE * fp = popen ("gnuplot 2>&1", "w");
  fprintf (fp, 
           "k(t)=a*exp(-b*t)*(1.-cos(c*t))\n"
           "a = 3e-4\n"
           "b = 1.5\n"
	   "\n"
           "D = %g\n"
           "n = 2.\n"
           "sigma = 1.\n"
           "rhol = 1.\n"
           "rhog = 1./1000.\n"
           "r0 = D/2.\n"
           "omega0 = sqrt((n**3-n)*sigma/((rhol+rhog)*r0**3))\n"
	   "\n"
           "c = 2.*omega0\n"
           "fit k(x) 'k-%d' via a,b,c\n"
	   "level = %d\n"
	   "res = D/%g*2.**level\n"
           "print \"fit \", res, a, b, c, D\n"
	   "\n"
	   "set table 'fit-%d'\n"
	   "plot [0:1] 2.*a*exp(-b*x)\n"
	   "unset table\n"
	   "\n"
	   "set print 'error' append\n"
	   "print res, c/2./omega0-1., D\n"
	   "\n"
	   "set print 'laplace' append\n"
	   "empirical_constant = 30.\n"
	   "print res, (1./(b**2.*D**3.))*empirical_constant**2, D\n"
	   "\n",
	   D, LEVEL, LEVEL, L0, LEVEL);
  pclose (fp);
}

#if QUADTREE
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D -s oscillation.gfv", "w");
  output_gfs (fp, t = t);
}

event adapt (i++) {
  adapt_wavelet ({c,u}, (double[]){5e-3,1e-3,1e-3}, LEVEL,
		 list = {p,u,pf,uf,g,c,alpha}, listb = {u,pf,uf,g,c,alpha});
}
#endif

/**
~~~gnuplot Evolution of the kinetic energy as a function of time for the spatial resolutions (number of grid points per diameter) indicated in the legend. The black lines are fitted decreasing exponential functions.
set output 'k.png'
set xlabel 'Time'
set ylabel 'Kinetic energy'
set logscale y
plot [0:1][8e-5:]'k-7' t "51.2" w l,				\
  'k-6' t "25.6" w l, 'k-5' t "12.8" w l,			\
  'fit-7' t "" w l lt 7, 'fit-6' t "" w l lt 7,			\
  'fit-5' t "" w l lt 7
~~~

~~~gnuplot Relative error in the oscillation frequency as a function of resolution.
set output 'frequency.png'
set xlabel 'Diameter (grid points)'
set ylabel 'Frequency error (%)'
set logscale x 2
unset grid
set xzeroaxis
set key spacing 1.5 top right
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'error' u (log($1)):(log(abs($2)*100.)) via a,b
plot 'error' u ($1):(abs($2)*100.) t "" w p pt 5 ps 2, \
  exp(f(log(x))) t ftitle(a,b)
~~~

~~~gnuplot Equivalent Laplace number estimated from the numerical damping of kinetic energy.
set output 'laplace.png'
set xlabel 'Diameter (grid points)'
set ylabel 'Equivalent Laplace number'
set grid
plot 'laplace' t "" w p pt 5 ps 2
~~~
*/
