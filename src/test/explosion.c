/**
# Two-dimensional explosion

We solve the Euler equations for a compressible gas. */

#include "compressible.h"

/**
We make boundary conditions free outflow. */

w.y[top]    = neumann(0);
w.y[bottom] = neumann(0);
w.x[left]   = neumann(0);
w.x[right]  = neumann(0);

/**
The domain spans $[-1:1]\times[-1:1]$. */

#define LEVEL 7

int main() {
  origin (-1, -1);
  size (2.);
  init_grid (1 << LEVEL);
  run(); 
}

/**
Initial conditions come from Toro's book (Riemann Solvers and
Numerical Methods for Fluid Dynamics, 3rd Edition, Springer Ed.)
Chapter 17 section 17.1.1 are given in terms of density ($\rho$),
pression ($p$), velocity ($u$) both at the left and right side of the
discontinuity placed at $R=0.4$. */

event init (t = 0)
{
  double R = 0.4 ;
  double rhoL = 1., rhoR = 0.125 ;
  double VmL = 0.0, VmR = 0.0 ;
  double pL = 1.0,  pR = 0.1 ;

  /**
  Left and right initial states for $\rho$, $\mathbf{w}$ and energy
  $E = \rho \mathbf{u}^2/2 + p/(\gamma-1)$. */

  foreach() {
    double r = sqrt(sq(x) + sq(y));
    double p;
    if (r <= R) {
      rho[] = rhoL;
      w.x[] = w.y[] = VmL;
      p = pL;
    }
    else {
      rho[] = rhoR;
      w.x[] = w.y[] = VmR;
      p = pR;
    }
    E[] = rho[]*sq(w.x[])/2. + p/(gammao - 1.);
    w.x[] *= x*rho[]/r;
    w.y[] *= y*rho[]/r;
  }
}

event print (t = 0.25)
{

  /**
  At $t=0.25$ we output the values of $\rho$ and the normal velocity
  $\mathbf{u}_n$ as functions of the radial coordinate. */

  foreach() {
    double r = sqrt(sq(x) + sq(y));
    double wn = (w.x[]*x + w.y[]*y)/r;
    printf ("%g %g %g\n", r, rho[], wn/rho[]);
  }

  /**
  For reference we also output a cross-section at $y=0$. */

  for (double x = 0; x <= 1; x += 1e-2)
    fprintf (stderr, "%g %.4f %.4f\n", x,
	     interpolate (rho, x, 0.),
	     interpolate (w.x, x, 0.));
}

/**
On quadtrees, we adapt the mesh by controlling the error on the density
field. */

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({rho}, (double[]){1e-5}, LEVEL + 1);
}
#endif

/**
## Results

Results are presented in terms of $\rho$ and normal velocity $u_n$ for
Cartesian (7 levels) and adaptive (8 levels) computations. The
numerical results compare very well with Toro's numerical experiments.

~~~gnuplot Radial density profile
set xrange [0:1]
set xlabel 'r'

set output 'rho.png'
set ylabel 'rho'
plot './out' u 1:2 w p pt 7 ps 0.2 t 'Adaptive', \
     './cout' u 1:2 w p pt 7 ps 0.2 t 'Cartesian'
~~~

~~~gnuplot Normal velocity
set output 'velocity.png'
set ylabel 'Normal velocity'
plot './out' u 1:3 w p pt 7 ps 0.2 t 'Adaptive', \
     './cout' u 1:3 w p pt 7 ps 0.2 t 'Cartesian'
~~~ 
*/
