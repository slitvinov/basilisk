/**
# Decaying two-dimensional turbulence

We solve the two-dimensional incompressible Euler equations using a
vorticity--streamfunction formulation. */

#include "grid/multigrid.h"
#include "navier-stokes/stream.h"

/**
The domain is square of size unity by default. The resolution is
constant at $256^2$. */

int main() {
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
$t=1000$. We fix the colorscale to $[-0.3:0.3]$.

![Evolution of the vorticity](turbulence/omega.gif) */

event output (i += 8; t <= 1000) {
  static FILE * fp = popen ("ppm2gif --delay 10 > omega.gif", "w");
  output_ppm (omega, fp, min = -0.3, max = 0.3);
}
