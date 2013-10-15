/**
# A streamfunction--vorticity solver for the Navier--Stokes equations

In two dimensions the incompressible, constant-density Navier--Stokes
equations can be written
$$
\partial_t\omega + \mathbf{u}\cdot\nabla\omega + \nu\nabla^2\omega = 0
$$
$$
\nabla^2\psi = \omega
$$
with $\nu$ the viscosity coefficient. The vorticity $\omega$ and
streamfunction $\psi$ are defined as
$$
\omega = \partial_x u_y - \partial_y u_x
$$
$$
u_x = - \partial_y\psi
$$
$$
u_y = \partial_x\psi
$$
The equation for the vorticity is an advection--diffusion equation
which can be solved using the flux--based advection scheme in
[`advection.h`](/src/advection.h). The equation for the streamfunction
is a [Poisson equation](/src/poisson.h). */

#include "advection.h"
#include "poisson.h"

/**
The field advected by the [advection solver](/src/advection.h) is
called `f`. We define $\omega$ as an alias for `f` for clarity. We
allocate the streamfunction field $\psi$ and a structure to store the
statistics on the convergence of the Poisson solver. */

#define omega f
scalar psi[];
mgstats mgpsi;

/**
Here we set the default boundary conditions for the
streamfunction. The default convention in Basilisk is no-flow through
the boundaries of the domain, i.e. they are a streamline
i.e. $\psi=$constant on the boundary. */

psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);

/**
We set the defaults values for the new field $\psi$ and for the `CFL`
(the default in `advection.h` is 0.5). This is done once at the
beginning of the simulation. */

event defaults (i = 0)
{
  CFL = 0.8;
  foreach()
    psi[] = 0.;
  boundary ({psi});
}

/**
At every timestep, but after all the other events for this timestep
have been processed (the '`last`' keyword), we update the streamfunction
field $\psi$ by solving a Poisson equation with the updated vorticity
field $\omega$ (which has just been advected/diffused). */

event streamfunction (i++, last)
{
  mgpsi = poisson (psi, omega);

/**
Using the new streamfunction, we can then update the components of the
velocity field. Since they are staggered relative to the
streamfunction, we need to average the horizontal and vertical
gradients, which gives the discrete expressions below. */

  trash ({u});
  foreach_face(x)
    u.x[] = (psi[0,-1] + psi[-1,-1] - psi[0,1] - psi[-1,1])/(4.*Delta);
  foreach_face(y)
    u.y[] = (psi[1,0] + psi[1,-1] - psi[-1,0] - psi[-1,-1])/(4.*Delta);
  boundary_normal ({u});
  boundary_tangent ({u});
}
