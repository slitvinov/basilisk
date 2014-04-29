/**
# Decaying two-dimensional turbulence

We solve the two-dimensional incompressible Euler equations using a
vorticity--streamfunction formulation. */

#include "grid/multigrid.h"
#include "navier-stokes/stream.h"

/**
The domain is square, of size unity and centered on the origin. The
resolution is constant at $256^2$. */

int main() {
  origin (-0.5, -0.5);
  init_grid (256);
  run();
}

/**
The initial condition for vorticity is just a white noise in the range
$[-1:1]$ .*/

event init (i = 0) {
  foreach()
    omega[] = noise();
}

/**
We generate images of the vorticity field every 8 timesteps up to
$t=1000$. We fix the colorscale to $[-0.3:0.3]$. */

event output (i += 8; t <= 1000) {
  output_ppm (omega, min = -0.3, max = 0.3);
}

/**
Finally, after completion of the simulation, we create a compressed
GIF animation. */

event gif (t = end) {

  /**
  We first convert each PPM image into GIF... */
  
  system ("convert out out-%04d.gif && rm -f out;"

	  /**
	  ... then use `gifsicle` to create the final compressed
	  animated GIF. */
	  
	  "gifsicle --colors 256 --optimize --delay 10"
	  " --loopcount=0 out-*.gif > omega.gif && rm -f out-*.gif");
}

/**
![Evolution of the vorticity](turbulence/omega.gif) */
