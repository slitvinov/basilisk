/**
# Rising bubble

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can solve either the
axisymmetric or planar version. */

#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif
#include "navier-stokes/centered.h"
// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"

#define LEVEL 8

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  size (2);
  init_grid (1 << LEVEL);
  
  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */

  rho1 = 1000., mu1 = 10.;
#if 1
  rho2 = 100., mu2 = 1., f.sigma = 24.5;
#else
  rho2 = 1., mu2 = 0.1, f.sigma = 1.96;
#endif

  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {

  /**
  The domain is a rectangle. We only simulate half the bubble. */
  
  mask (y > 0.5 ? top : none);

  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */

  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume. */

event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
  }
  printf ("%g %g %g %g %g %g %d %d %d\n", 
	  t, sb, -1., xb/sb, vb/sb, dt, mgp.i, mgpf.i, mgu.i);
}

/**
At $t=3$ we output the shape of the bubble. */

event interface (t = 3.) {
  output_facets (f, stderr);
}

/**
If gfsview is installed on the system, we can also visualise the
simulation as it proceeds. */

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview2D rising.gfv", "w");
  output_gfs (fp);
}
#endif

/**
## Results

The final shape of the bubble is compared to that obtained with the
MooNMD Lagrangian solver (see [the FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html))
at the highest resolution. We also display the shape of the
axisymmetric version of the test. The axisymmetric bubble moves much
faster.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set size ratio -1
set grid
plot [][0:0.4]'../c1g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              'log' u 1:2 w l t 'Basilisk', \
              '../rising-axi/log' u 1:2 w l t 'Basilisk (axisymmetric)'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time.
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c1g3l4.txt' u 1:5 w l t 'MooNMD', \
                 'out' u 1:5 w l t 'Basilisk', \
                 '../rising-axi/out' u 1:5 w l t 'Basilisk (axisymmetric)'
~~~

*/
