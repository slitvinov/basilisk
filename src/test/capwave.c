/**
# Capillary wave

This is the classical test case first proposed in [Popinet & Zaleski,
1999](/src/references.bib#popinet1999).

We use a constant-resolution grid, the Navier--Stokes solver with VOF
interface tracking and surface tension. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

/**
The interface is represented by the volume fraction field *c*. There
are no diffusive tracers. */

scalar c[], * interfaces = {c}, * tracers = NULL;

/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.x[left]   = dirichlet(0);
uf.x[right]  = dirichlet(0);
uf.y[top]    = dirichlet(0);
uf.y[bottom] = dirichlet(0);

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */

double se = 0, ne = 0;

int main() {

  /**
  The domain is 2x2 to minimise finite-size effects. The surface
  tension is one and the viscosity is constant. */

  L0 = 2.;
  Y0 = -L0/2.;
  c.sigma = 1.;
  TOLERANCE = 1e-6;
  const face vector muc[] = {0.0182571749236, 0.0182571749236};
  mu = muc;

  /**
  We vary the resolution to check for convergence. */

  for (N = 16; N <= 128; N *= 2) {
    se = ne = 0;
    run();
  }
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - 0.01*cos (2.*pi*x));
  fractions (phi, c);
}

/**
We output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 3.04290519077e-3; t <= 2.2426211256) {

  /**
  To get an accurate amplitude, we reconstruct the height function
  field and take the corresponding maximum. */

  vector h[];
  heights (c, h);
  double max = - HUGE;;
  foreach() 
    if (c[] > 0 && c[] < 1) {
      double yi = y + h.y[];
      if (yi > max)
	max = yi;
    }

  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points *N*. */

  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*11.1366559937, max);
  fflush (fp);

  /**
  To compute the RMS error, we get data from the reference file
  *prosperetti* and add the difference to the accumulated error. */

  static FILE * fp1 = fopen ("../prosperetti", "r");
  double t1, max1;
  fscanf (fp1, "%lf %lf", &t1, &max1);
  se += sq(max - max1); ne++;
}

/**
At the end of the simulation, we output on standard error the
resolution (number of grid points per wavelength) and the relative RMS
error. */

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0, sqrt(se/ne)/0.01);

#if QUADTREE
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D -s oscillation.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
## Results

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

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html)
*/
