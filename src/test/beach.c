/**
# Solitary wave run-up on a plane beach

We use the [Green-Naghdi](/src/green-naghdi.h) solver to reproduce
this classical test case based on the experiments of [Synolakis,
1987](/src/references.bib#synolakis1987). */

#include "grid/multigrid1D.h"
#include "green-naghdi.h"

/**
The problem is non-dimensionalised by the water depth $h_0$ and the
acceleration of gravity. The amplitude of the solitary wave is 0.28
and the slope of the beach is 1/19.85. The length $L$ is the estimated
wavelength of the solitary wave (see section 4.3 of [Yamazaki et al,
2009](/src/references.bib#yamazaki2009)). We set the coordinate system
(*X0* and *L0*) so that the origin is the intersection of the beach
with the water level. */

double h0 = 1., a = 0.28, L;
double slope = 1./19.85;

int main() {
  double k = sqrt(3.*a/4/cube(h0));
  L = 2./k*acosh(sqrt(1./0.05));
  X0 = - h0/slope - L/2. - L;
  L0 = 6.*L;
  N = 1024;
  G = 1.;
  alpha_d = 1.;
  
  /**
  We try to tune the "breaking slope" to get better agreement for
  $t=20$. */
  
  breaking = 0.4;
  run();
}

/**
The initial wave is the [analytical soliton](soliton.c) for the
Green-Naghdi equations. */

double sech2 (double x) {
  double a = 2./(exp(x) + exp(-x));
  return a*a;
}

double soliton (double x, double t)
{
  double c = sqrt(G*(1. + a)*h0), psi = x - c*t;
  double k = sqrt(3.*a*h0)/(2.*h0*sqrt(h0*(1. + a)));
  return a*h0*sech2 (k*psi);
}

/**
See figure 5 of [Yamazaki et al,
2009](/src/references.bib#yamazaki2009) for the definition of the
bathymetry. */

event init (i = 0)
{
  double c = sqrt(G*(1. + a)*h0);
  foreach() {
    double eta = soliton (x + h0/slope + L/2., t);
    zb[] = max (slope*x, -h0);
    h[] = max (0., eta - zb[]);
    u.x[] = c*eta/(h0 + eta);
  }
}

/**
Friction is important for this test case. We implement a simple
time-implicit quadratic bottom friction and tune the coefficient to
obtain a runup comparable with the experiment. */

event friction (i++) {
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 5e-3*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
}

/**
We use gnuplot to display an animation while the simulation is
running. */

event gnuplot (i += 5) {
  static FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p [-20:12][-0.2:0.6]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 t '' w filledcu lc -1\n", t);
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fprintf (stderr, "%.3f %.3f\n", t, statsf(u.x).max);
}

/**
We output profiles at the same times as the experimental data, in
separate files indexed by the time. */

event output (t <= 65; t += 5) {
  char name[80];
  sprintf (name, "out-%g", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "\n");
  fclose (fp);
}

/**
The agreement with the experiments (circles) is satisfactory and is
comparable to the numerical results of [Yamazaki et al,
2009](/src/references.bib#yamazaki2009) Figure 6, obtained with a
different set of depth-averaged equations and the results of [Bonneton
et al, 2011](/src/references.bib#bonneton2011), Figure 8, although
this latter model seems to do a better job of capturing breaking
around $t=20$ (they use a more sophisticated breaking criterion).

~~~gnuplot Comparison of model predictions and experimental snapshots
set term @PNG enhanced size 640,640 font ",8"
set xrange [-20:12]
set yrange [-0.2:0.6]
set ytics -0.2,0.2,0.6
set key top left
set multiplot layout 6,2 scale 1.05,1.1
unset xtics
plot 'out-10' w l lw 2 t 't=10', '' u 1:3 w filledcu lc -1 t '', \
     '../t-10' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-15' w l lw 2 t 't=15', '' u 1:3 w filledcu lc -1 t '', \
     '../t-15' pt 6 lc -1 ps 0.5 t ''
set ytics -0.2,0.2,0.6
plot 'out-20' w l lw 2 t 't=20', '' u 1:3 w filledcu lc -1 t '', \
     '../t-20' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-25' w l lw 2 t 't=25', '' u 1:3 w filledcu lc -1 t '', \
     '../t-25' pt 6 lc -1 ps 0.5 t ''
set ytics -0.2,0.2,0.6
plot 'out-30' w l lw 2 t 't=30', '' u 1:3 w filledcu lc -1 t '', \
     '../t-30' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-35' w l lw 2 t 't=35', '' u 1:3 w filledcu lc -1 t '', \
     '../t-35' pt 6 lc -1 ps 0.5 t ''
set ytics -0.2,0.2,0.6
plot 'out-40' w l lw 2 t 't=40', '' u 1:3 w filledcu lc -1 t '', \
     '../t-40' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-45' w l lw 2 t 't=45', '' u 1:3 w filledcu lc -1 t '', \
     '../t-45' pt 6 lc -1 ps 0.5 t ''
set ytics -0.2,0.2,0.6
plot 'out-50' w l lw 2 t 't=50', '' u 1:3 w filledcu lc -1 t '', \
     '../t-50' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-55' w l lw 2 t 't=55', '' u 1:3 w filledcu lc -1 t '', \
     '../t-55' pt 6 lc -1 ps 0.5 t ''
set ytics -0.2,0.2,0.6
set xtics
plot 'out-60' w l lw 2 t 't=60', '' u 1:3 w filledcu lc -1 t '', \
     '../t-60' pt 6 lc -1 ps 0.5 t ''
unset ytics
plot 'out-65' w l lw 2 t 't=65', '' u 1:3 w filledcu lc -1 t '', \
     '../t-65' pt 6 lc -1 ps 0.5 t ''
unset multiplot
~~~
*/
