/**
# Incompressible Navier--Stokes solver (MAC formulation)

We wish to approximate numerically the incompressible Navier--Stokes
equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left(-\nabla p + \nabla\cdot(\mu\nabla\mathbf{u})\right)
$$
$$
\nabla\cdot\mathbf{u} = 0
$$

We will use the generic time loop, a CFL-limited timestep and we will
need to solve a Poisson problem. */

#include "run.h"
#include "timestep.h"
#include "poisson.h"

/**
The Markers-And-Cells (MAC) formulation was first described in the
pioneering paper of [Harlow and Welch,
1964](/src/references.bib#harlow1964). It relies on a *staggered*
discretisation of the velocity components `u.x` and `u.y`, relative to
the (centered) pressure `p`. This guarantees the consistency of the
discrete gradient, divergence and Laplacian operators and leads to a
stable (mode-free) integration. */

scalar p[];
staggered vector u[];

/**
The parameters are the viscosity coefficient $\mu$ and the specific
volume $\alpha = 1/\rho$ (with default unity). $\alpha$ is a staggered
vector because we will need its values at the staggered locations of
the velocity components. 

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

double mu = 0.;
staggered vector alpha;
mgstats mgp;

/**
The default boundary conditions are slip walls. */

u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

/**
The default velocity and pressure are zero. */

event defaults (i = 0)
{
  foreach_face()
    u.x[] = 0.;
  foreach()
    p[] = 0.;
  boundary ({p,u});
}

/**
We apply boundary conditions again after user initialisation. */

event init (i = 0)
{
  boundary ({p,u});
}

/**
## Time integration

### Advection--Diffusion 

In a first step, we compute $\mathbf{u}_*$
$$
\frac{\mathbf{u}_* - \mathbf{u}_n}{dt} = \nabla\cdot\mathbf{S}
$$
with $\mathbf{S}$ the symmetric tensor
$$
\mathbf{S} = - \mathbf{u}\otimes\mathbf{u} + \mu\nabla\mathbf{u} =
\left(\begin{array}{cc}
- u_x^2 + 2\mu\partial_xu_x & - u_xu_y + \mu(\partial_yu_x + \partial_xu_y)\\
\ldots & - u_y^2 + 2\mu\partial_yu_y
\end{array}\right)
$$ */

event advance (i++)
{

/**
We allocate a local symmetric tensor field. To be able to compute the
divergence of the tensor at the staggered locations, we need to
compute the diagonal components at the center of cells and the
off-diagonal component at the vertices. */

  symmetric tensor S[];
  S.x.y.d.x = S.x.y.d.y = -1; // staggering
  trash ({S});

/**
We average the velocity components at the center to compute the
diagonal components. */

  foreach()
    foreach_dimension()
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*mu*(u.x[1,0] - u.x[])/Delta;

/**
We average horizontally and vertically to compute the off-diagonal
component at the vertices. */

  foreach_vertex()
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      mu*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/Delta;

/**
We only need to apply boundary conditions to the diagonal components. */

  boundary ({S.x.x, S.y.y});

/**
Finally we compute
$$
\mathbf{u}_* = \mathbf{u}_n + dt\nabla\cdot\mathbf{S}
$$ */

  foreach_face()
    u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/Delta;
}

/**
### Projection 

In a second step we compute
$$
\mathbf{u}_{n+1} = \mathbf{u}_* - \alpha\nabla p
$$
with the condition
$$
\nabla\cdot\mathbf{u}_{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\alpha\nabla p) = \nabla\cdot\mathbf{u}_*
$$ */

event projection (i++)
{
  
/**
We allocate a local scalar field and compute the divergence of
$\mathbf{u}_*$. */

  scalar div[];
  boundary ((scalar *){u});
  foreach()
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/Delta;

/**
We solve the Poisson problem. */

  mgp = poisson (p, div, alpha);

/**
And compute $\mathbf{u}_{n+1}$ using $\mathbf{u}_*$ and $p$. If
$\alpha$ is not defined we set it to one. */

  if (alpha.x)
    foreach_face()
      u.x[] -= alpha.x[]*(p[] - p[-1,0])/Delta;
  else
    foreach_face()
      u.x[] -= (p[] - p[-1,0])/Delta;
  boundary ((scalar *){u});

/**
Finally we obtain the timestep for the next iteration by applying the
CFL condition to the new velocity field. */

  dt = timestep (u);
}
