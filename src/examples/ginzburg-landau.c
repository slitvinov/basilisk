/**
# The complex Ginzburg--Landau equation

The complex [Ginzburg--Landau
equation](http://codeinthehole.com/static/tutorial/index.html)
$$
\partial_t A = A + \left( 1 + i \alpha \right) \nabla^2 A - \left( 1 + i
\beta \right)  \left| A \right|^2 A
$$
with $A$ a complex number, is a classical model for phenomena
exhibiting [Hopf
bifurcations](http://en.wikipedia.org/wiki/Hopf_bifurcation) such as
Rayleigh-Bénard convection or superconductivity.

Posing $A_r=Re(A)$ and $A_i=Im(A)$ one gets the
coupled reaction--diffusion equations.
$$
\partial_t A_r = \nabla^2 A_r + A_r  \left( 1 - \left| A \right|^2 \right)
   - \alpha \nabla^2 A_i + \left| A \right|^2 \beta A_i
$$
$$
\partial_t A_i = \nabla^2 A_i + A_i  \left( 1 - \left| A \right|^2 \right)
   + \alpha \nabla^2 A_r - \left| A \right|^2 \beta A_r
$$

This system can be solved with the reaction--diffusion solver. */

#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"

scalar Ar[], Ai[], A2[];

/**
In this example, we only consider the case when $\alpha=0$. */

double beta;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
mgstats mgd1, mgd2;

/**
## Parameters

We change the size of the domain `L0`. */

int main() {
  beta = 1.5;
  size (100);
  init_grid (256);
  run();
}

/**
## Initial conditions 

We use a white noise in $[-10^{-4}:10^{-4}]$ for both components. */

event init (i = 0) {
  foreach() {
    Ar[] = 1e-4*noise();
    Ai[] = 1e-4*noise();
  }
  boundary ({Ar,Ai});
}

/**
## Time integration */

event integration (i++) {

  /**
  We first set the timestep according to the timing of upcoming
  events. We choose a maximum timestep of 0.05 which ensures the stability
  of the reactive terms for this example. */

  dt = dtnext (t, 0.05);

  /**
  We compute $|A|^2$. */

  foreach()
    A2[] = sq(Ar[]) + sq(Ai[]);
  boundary ({A2});

  /**
  We use the diffusion solver (twice) to advance the system from $t$
  to $t+dt$. */

  scalar r[], lambda[];
  foreach() {
    r[] = A2[]*beta*Ai[];
    lambda[] = 1. - A2[];
  }
  mgd1 = diffusion (Ar, dt, r = r, beta = lambda);
  foreach() {
    r[] = - A2[]*beta*Ar[];
    lambda[] = 1. - A2[];
  }
  mgd1 = diffusion (Ai, dt, r = r, beta = lambda);
}

/**
## Outputs

Here we create mpeg animations for both components. The `spread`
parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movies (i += 3; t <= 150) {
  fprintf (stderr, "%g %g\n", t, sqrt(normf(A2).max));

  static FILE * fp = popen ("ppm2mpeg > Ai.mpg", "w");
  output_ppm (Ai, fp, spread = 2, linear = true);
  static FILE * fp1 = popen ("ppm2mpeg > A2.mpg", "w");
  output_ppm (A2, fp1, spread = 2, linear = true);
}

/**
At the end of the simulation, we create snapshot images of both
fields, in PNG format. */

event pictures (t = end) {
  output_ppm (Ai, file = "Ai.png", spread = 2, linear = true);
  output_ppm (A2, file = "A2.png", spread = 2, linear = true);
}

/**
For the value of $\beta$ we chose, we get a typical "frozen state"
composed of "cellular structures" for $|A|^2$ and stationary spirals
for $A_i$.

|:------:|:-----:|
| [![](ginzburg-landau/A2.png)](ginzburg-landau/A2.mpg) | [![](ginzburg-landau/Ai.png)](ginzburg-landau/Ai.mpg) |
| $|A|^2$                     |    $A_i$                    |

  : Final states (click on images for animations)

## See also

* [The Brusselator](brusselator.c) */
