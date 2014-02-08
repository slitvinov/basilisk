/**
# Merging of two vertices

Studying the interaction of two incompressible vortices is interesting
for example as a basic model of two-dimensional turbulence. Here we
solve the incompressible 2D Euler equations using a
vorticity--streamfunction formulation. */

#include "navier-stokes/stream.h"

/**
The domain is centered on $(0,0)$ and the maximum level of refinement
is 8 i.e. the initial grid has $N=2^8=256$ grid points per
dimension. */

#define MAXLEVEL 8

int main()
{
  X0 = Y0 = -0.5;
  N = 1 << MAXLEVEL;
  run();
}

/**
The initial vorticity field is composed of two Gaussian vortices
separated by twice *dd* and with caracteristic radii *dd/10*. */

event init (i = 0)
{
  double dd = 0.1;
  foreach()
    omega[] = (exp(-(sq(x - dd) + sq(y))/(dd/10.)) +
	       exp(-(sq(x + dd) + sq(y))/(dd/10.)));
}

/**
We output some statistics on the vorticity field and Poisson solver at
the start and end of the simulation. */

event logfile (t = {0,30}) {
  stats s = statsf (omega);
  fprintf (stderr, "%g %d %g %g %d\n", t, i, dt, s.sum, mgpsi.i);
}

/**
We output the vorticity and level fields at regular intervals in a
format compatible with gnuplot. */

event output (t += 5) {
  static int nf = 0;
  printf ("file: omega-%d\n", nf);
  output_field ({omega}, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({l});
  nf++;
}

/**
If we are using a quadtree grid, it is adapted using wavelet error
control on $\omega$. */

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({omega}, (double[]){1e-2}, MAXLEVEL, list = {omega, psi});
}
#endif

/**
## Results

After running and processing by gnuplot (using [stream.plot]()) we get:

![Evolution of the vorticity field with time.](stream/plot.png)

![Evolution of level of refinement with time.](stream/level.png)

## See also

* [Merging of two vortices (centered Euler solver)](vortex.c)
* [Coalescence of a pair of Gaussian vortices (Gerris logo)](http://gerris.dalembert.upmc.fr/gerris/examples/examples/logo.html) */
