/**
# Rising bubble

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. */

#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

/**
Hysing et al. consider two cases (1 and 2), with the densities, dynamic
viscosities and surface tension of fluid 1 and 2 given below. */

#define rho1 1000.
#define mu1 10.

#if 1
# define rho2 100.
# define mu2 1.
# define SIGMA 24.5
#else
# define rho2 1.
# define mu2 0.1
# define SIGMA 1.96
#endif

#define LEVEL 8

/**
The interface between the two-phases will be tracked with the volume
fraction field *f*. We allocate the fields for the variable density
and viscosity (*alphav*, *alphacv* and *muv* respectively). */

scalar f[];
scalar * interfaces = {f};
face vector alphav[];
scalar alphacv[];
face vector muv[];

/**
The boundary conditions are slip lateral walls and no-slip on the top
and bottom walls. */

u.t[top]     = dirichlet(0);
uf.t[top]    = dirichlet(0);
u.t[bottom]  = dirichlet(0);
uf.t[bottom] = dirichlet(0);

uf.n[left]  = 0;
uf.n[right] = 0;

/**
The domain will span $[0.5:1]\times[0:2]$ and will be resolved with
$64\times 256$ grid points. */

int main() {
  size (2);
  init_grid (1 << LEVEL);
  run();
}

event init (t = 0) {
  
  /**
  The density and viscosity are defined by the variable fields we
  allocated above. We also set the surface tension for interface
  *f*. */
  
  alpha = alphav;
  alphac = alphacv;
  mu = muv;
  f.sigma = SIGMA;

  /**
  The domain is a rectangle. We only simulate half the bubble. */
  
  mask (x > 1. ? right :
	x < 0.5 ? left :
	none);

  /**
  The bubble is centered on (0.5,0.5) and has a radius of 0.25. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x - 0.5) + sq(y - 0.5) - sq(0.25);
  fractions (phi, f);
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 0.98;
}

/**
The density and viscosity are defined using the arithmetic average. */

#define rho(f) ((f)*rho1 + (1. - (f))*rho2)
#define mu(f)  ((f)*mu1 + (1. - (f))*mu2)

event properties (i++) {
  foreach_face() {
    double fm = (f[] + f[-1,0])/2.;
    alphav.x[] = 1./rho(fm);
    muv.x[] = mu(fm);
  }
  foreach()
    alphacv[] = 1./rho(f[]);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume. */

event logfile (i++) {
  double yb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:yb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    yb += y*dv;
    vb += u.y[]*dv;
    sb += dv;
  }
  printf ("%g %g %g %g %g %g %d %d %d\n", 
	  t, sb, -1., yb/sb, vb/sb, dt, mgp.i, mgpf.i, mgu.i);
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
  static FILE * fp = popen("gfsview2D ../rising.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
## Results

The final shape of the bubble is compared to that obtained with the
MooNMD Lagrangian solver (see [the FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html))
at the highest resolution.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set size ratio -1
set grid
plot [0.1:0.9][0.8:1.4]'../c1g3l4s.txt' w l t 'MooNMD', 'log' w l t 'Basilisk'
~~~

The agreement for the bubble raise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time.
reset
set grid
set xlabel 'Time'
plot [0:3][0:0.3]'../c1g3l4.txt' u 1:5 w l t 'MooNMD', \
                 'out' u 1:5 w l t 'Basilisk
~~~

*/
