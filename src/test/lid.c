/**
# Lid-driven cavity at Re=1000

We use the multigrid implementation (rather than the default quadtree
implementation) and the MAC Navier--Stokes solver. */

#include "grid/multigrid.h"
#include "navier-stokes1.h"

/**
Here we define the domain geometry: a square box of size unity
centered on (0,0). We also set the viscosity and some parameters 
controlling the numerical scheme. */

void parameters()
{ 
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // number of grid points
  N = 64;
  // viscosity
  NU = 1e-3;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

/**
The default boundary conditions are symmetry (i.e. slip walls). We
need no-slip on three boundaries and $u=1$ on the top
boundary. Boundary conditions are set by defining the value of the
'ghost cell' e.g. for the top boundary we have

![](/boundary.png)

with
$$
u_b = (u + u_g)/2
$$
this gives a ghost cell value for the top boundary of */

u.x[top]    = 2. - u.x[];

/**
For the other no-slip boundaries this gives */

u.x[bottom] = - u.x[];
u.y[left]   = - u.y[];
u.y[right]  = - u.y[];

/**
The default velocity is zero everywhere. We still need to define an
init() function. */

void init() {}

/**
We define an auxilliary function which computes the total kinetic
energy. */

static double energy()
{
  double se = 0.;
  foreach(reduction(+:se))
    se += (sq(u.x[] + u.x[1,0)] + sq(u.y[] + u.y[0,1))]/8.*delta*delta;
  return se;
}

/**
We want the simulation to stop when we are close to steady state. To
do this we store the `u.x` field of the previous timestep in an
auxilliary variable `un`. */

scalar un[];

event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

/**
Every 10 timesteps we check whether $u$ has changed by more than
10^-4^. If it has not, the event returns 1 which stops the
simulation. We also output the evolution of the kinetic energy on
standard error. */

event logfile (i += 10; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-4)
    return 1; /* stop */
  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

/**
Every 100 timesteps we output a binary representation of `u.x`
bilinearly-interpolated on an N x N grid. */

event outputfile (i += 100) output_matrix (u.x, stdout, N, linear = true);

/**
This function is called after completion of the simulation. We write
in the `xprof` and `yprof` files the interpolated values of `u.x` and
`u.y` along the two profiles. */

void end()
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}

/**
And finally we call the `run()` method of the Navier--Stokes solver. */

int main() { run (); }

/**
## Results

After processing by gnuplot (i.e. using `make lid/plot.png` with
[lid.plot]()) we get

![Horizontal profile of the $y$-component of the velocity on the
 centerline of the box.](lid/plot.png)

The results can be made closer to [Ghia et al, 1982](/src/references.bib#ghia82) 
by increasing the resolution. */
