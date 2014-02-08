/**
# Merging of two vortices (centered Euler solver)

This test is similar to [stream.c]() but uses the centered
Navier--Stokes solver (without viscosity). It also shows how to
convert a vorticity field into a velocity field. */

#include "navier-stokes/centered.h"

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
For the centered Navier--Stokes solver, the primary variables are the
velocity and pressure field. We need to convert the initial vorticity
field into the velocity field. To do so we first declare the
(temporary) streamfunction $\psi$ and vorticity $\omega$ fields. We
also set appropriate boundary conditions for the streamfunction. */

event init (t = 0)
{
  scalar psi[], omega[];

  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

/**
We then initialise both fields, using the same initial condition for
the vorticity as in [stream.c](), and apply boundary conditions. Note
that it is necessary to initialise the streamfunction, as the Poisson
solver requires an initial guess. */

  double dd = 0.1;
  foreach() {
    omega[] = (exp(-(sq(x - dd) + sq(y))/(dd/10.)) +
	       exp(-(sq(x + dd) + sq(y))/(dd/10.)));
    psi[] = 0.;
  }
  boundary ({psi,omega});

/**
We then solve the Poisson equation
$$
\nabla^2\psi = \omega
$$
and compute the centered velocity components by differentation of the
streamfunction i.e.
$$
u_x = - \partial_y\psi
$$
$$
u_y = \partial_x\psi
$$ */

  poisson (psi, omega);
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});
}

/**
For convenience, we also define a function which, given a velocity
field $\mathbf{u}$, fills a scalar field $\omega$ with the vorticity
field
$$
\omega = \partial_x u_y - \partial_y u_x
$$ */

void vorticity (vector u, scalar omega)
{
  foreach()
    omega[] = (u.y[1,0] - u.y[-1,0] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
  boundary ({omega});
}

/**
We output some statistics on the vorticity field and Poisson solver at
the start and end of the simulation. */

event logfile (t = {0,30}) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);
  fprintf (stderr, "%g %d %g %g %d\n", t, i, dt, s.sum, mgp.i);
}

/**
We make animations of the vorticity and level of refinement. */

event movie (t += 0.2; t <= 30) {
  static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, fp, linear = true);

  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  foreach()
    omega[] = level;
  output_ppm (omega, fp1, spread = 2);
}

/**
We output the vorticity and level fields at regular intervals in a
format compatible with gnuplot. */

event output (t += 5) {
  static int nf = 0;
  scalar omega[];
  vorticity (u, omega);
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
control on both components of the velocity field. Note that the error
thresholds need to be specified twice (once for each component of
vector $\mathbf{u}$). */

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-4,1e-4}, MAXLEVEL);
}
#endif

/**
## Results

After running and processing by gnuplot (using [vortex.plot]()) we get
the following pictures and animations.

![[Evolution of the vorticity field with time.](vortex/vort.mpg)](vortex/plot.png)

![[Evolution of level of refinement with time.](vortex/level.mpg)](vortex/level.png) */
