/**
# Shock reflection by a circular cylinder

The evolution of an initial "step" wave is modelled using the
Saint-Venant equations. The wave interacts with a circular cylinder
described using embedded solid boundaries. Adaptivity is used to track
the wave fronts. This example is discussed in [An and Yu,
2012](/src/references.bib#an2012). */

#include "saint-venant.h"

int LEVEL = 9;

/**
We define a new boundary for the cylinder. */

bid cylinder;

int main() {
  size (5.);
  G = 9.81;
  origin (-L0/2., -L0/2.);
  init_grid (1 << LEVEL);
  run();
}

/**
We impose height and velocity on the left boundary. */

#define H0 3.505271526
#define U0 6.29033769408481

h[left]   = H0;
eta[left] = H0;
u.n[left] = U0;

event init (i = 0) {

  /**
  The geometry is defined by masking and the initial step function is
  imposed. */
  
  mask (sq(x + 0.5) + sq(y) < sq(0.5) ? cylinder : none);
  foreach() {
    h[] = (x <= -1 ? H0 : 1.);
    u.x[] = (x <= -1 ? U0 : 0.);
  }
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

/**
We generate movies of depth and level of refinement. */

event movie (t += 0.0025) {
  static FILE * fp = popen ("ppm2mpeg > depth.mpg", "w");
  output_ppm (h, fp,
	      min = 0.1, max = 6, map = cool_warm, n = 400, linear = true);
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, map = cool_warm, min = 4, max = LEVEL, n = 400);
}

/**
At $t=0.3$ we generate snapshots of the depth and level of refinement... */

event image (t = 0.3) {
  output_ppm (h, file = "depth.png",
	      min = 0.1, max = 6, map = cool_warm, n = 400, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "level.png", map = cool_warm,
	      min = 4, max = LEVEL, n = 400);
}

/**
...which gives the following results

![[Animation](shock/depth.mpg) of the fluid depth. Snapshot at
 $t=0.3$.](shock/depth.png)

![[Animation](shock/level.mpg) of the level of refinement. Snapshot at
 $t=0.3$.](shock/level.png)

The mesh is adapted according to the error on the height field. */

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-2}, LEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

/**
## See also

* [Same case with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/shock.html)
*/
