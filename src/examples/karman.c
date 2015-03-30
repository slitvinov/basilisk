/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

We use the centered Navier-Stokes solver and advect the passive tracer
*f*. */

#include "navier-stokes/centered.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  run(); 
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
We add a new boundary condition for the cylinder. The tangential
velocity on the cylinder is set to zero. */

bid cylinder;
u.t[cylinder] = dirichlet(0.);

event init (t = 0) {

  /**
  To make a long channel, we set the *top* boundary for $y > 0.5$ and
  the *bottom* boundary for $y < -0.5$. The *cylinder* has a radius of
  0.0625. */

  mask (y >  0.5 ? top :
	y < -0.5 ? bottom :
	sq(x) + sq(y) < sq(0.0625) ? cylinder :
	-1);

  /**
  We set a constant viscosity corresponding to a Reynolds number of
  160, based on the cylinder diameter (0.125) and the inflow velocity
  (1). We also set the initial velocity field and tracer
  concentration. */
  
  const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  foreach() {
    u.x[] = 1.;
    f[] = 0.;
  }
  boundary ({f});
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 2; t <= 15.) {
  static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
  scalar vorticity[];
  foreach()
    vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
  boundary ({vorticity});
  output_ppm (vorticity, fp, box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true);
  static FILE * fp1 = popen ("ppm2mpeg > f.mpg", "w");
  output_ppm (f, fp1, box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = true, min = 0, max = 1);
}

/**
... and the corresponding snapshots at the end of the simulation. This
gives the following figures and animations.

![[Animation](karman/vort.mpg) of the vorticity field.](karman/vort.png)

![[Animation](karman/f.mpg) of the tracer field.](karman/f.png) 
*/

event snapshots (t = end)
{
  scalar vorticity[];
  foreach()
    vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
  boundary ({vorticity});
  output_ppm (vorticity, file = "vort.png",
	      box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true);
  output_ppm (f, file = "f.png", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = true, min = 0, max = 1);
}

/**
If gfsview is installed on your system you can use this to visualise
the simulation as it runs. */

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s ../karman.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
We adapt according to the error on the velocity and tracer fields. */

event adapt (i++) {
  adapt_wavelet ({u,f}, (double[]){3e-2,3e-2,3e-2}, 9, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
