/**
# Incompressible Navier--Stokes solver (centered formulation)

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
#include "bcg.h"
#include "viscosity.h"

/**
The Markers-And-Cells (MAC) formulation was first described in the
pioneering paper of [Harlow and Welch,
1965](/src/references.bib#harlow1965). It relies on a *face*
discretisation of the velocity components `u.x` and `u.y`, relative to
the (centered) pressure `p`. This guarantees the consistency of the
discrete gradient, divergence and Laplacian operators and leads to a
stable (mode-free) integration. */

scalar p[];
vector u[];
scalar pf[];
face vector uf[];

/**
The parameters are the viscosity coefficient $\mu$ and the specific
volume $\alpha = 1/\rho$ (with default unity). $\alpha$ is a face
vector because we will need its values at the face locations of
the velocity components. 

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

face vector alpha; // default unity
face vector mu;    // default zero
mgstats mgp, mgpf, mgu;

/**
The default velocity and pressure are zero. */

// fixme: (normal) BCs for uf needs to be consistent with (normal) BCs for u

event defaults (i = 0)
{
  const face vector unity[] = {1,1};
  alpha = unity;

  CFL = 0.8;
  foreach() {
    foreach_dimension()
      u.x[] = 0.;
    p[] = pf[] = 0.;
  }
  boundary ({p,pf,u});
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

void prediction()
{
  scalar gx[], gy[];
  vector g;
  g.x = gx; g.y = gy;
  foreach() {
    if (u.x.gradient)
      foreach_dimension()
	g.x[] = u.x.gradient (u.x[-1,0], u.x[], u.x[1,0])/Delta;
    else
      foreach_dimension()
	g.x[] = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
  }
  boundary ((scalar *){g});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1,0])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*Delta/2.;
    double fyy = u.y[i,0] < 0. ? u.x[i,1] - u.x[i,0] : u.x[i,0] - u.x[i,-1];
    uf.x[] -= dt*u.y[i,0]*fyy/(2.*Delta);
  }
  boundary ((scalar *){uf});
}

// fixme: share with MAC solver
mgstats project (face vector u, scalar p, (const) face vector alpha)
{
  scalar div[];
  foreach()
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/Delta;
  mgstats mgp = poisson (p, div, alpha);
  foreach_face()
    u.x[] -= alpha.x[]*(p[] - p[-1,0])/Delta;
  boundary_normal ({u});
  boundary_tangent ({u});
  return mgp;
}

event stability (i++) {}

event advance (i++)
{
  prediction();
  mgpf = project (uf, pf, alpha);

  (const) face vector alphaf = alpha;
  foreach_dimension() {
    vector flux[];
    tracer_fluxes (u.x, uf, flux, dt);
    foreach()
      u.x[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/Delta +
      dt*(uf.x[1,0]*alphaf.x[1,0]*(p[1,0] - p[]) +
	  uf.y[0,1]*alphaf.y[0,1]*(p[1,0] - p[-1,0] + p[1,1] - p[-1,1])/2. -
	  uf.x[]*alphaf.x[]*(p[] - p[-1,0]) -
	  uf.y[]*alphaf.y[]*(p[1,0] - p[-1,0] + p[1,-1] - p[-1,-1])/2.)
      /(2.*sq(Delta));
  }

  if (mu.x && dt > 0.) {
    scalar alphac[];
    foreach() {
      alphac[] = 0.;
      foreach_dimension() {
        u.x[] -= (alphaf.x[]*(p[] - p[-1,0]) +
		  alphaf.x[1,0]*(p[1,0] - p[]))/(2.*Delta);
	alphac[] += alphaf.x[] + alphaf.x[1,0];
      }
      alphac[] *= dt/4;
    }
    boundary ((scalar *){u});
    mgu = viscosity (u, mu, alphac);
    foreach()
      foreach_dimension()
        u.x[] += (alphaf.x[]*(p[] - p[-1,0]) +
		  alphaf.x[1,0]*(p[1,0] - p[]))/(2.*Delta);
  }
  
  boundary ((scalar *){u});
}

event projection (i++)
{
  trash ({uf});
  foreach_face()
    uf.x[] = (u.x[] + u.x[-1,0])/2.;
  boundary_normal ({uf});
  mgp = project (uf, p, alpha);

  (const) face vector alphaf = alpha;
  foreach()
    foreach_dimension()
      u.x[] -= (alphaf.x[]*(p[] - p[-1,0]) +
		alphaf.x[1,0]*(p[1,0] - p[]))/(2.*Delta);
  boundary ((scalar *){u});

/**
Finally we obtain the timestep for the next iteration by applying the
CFL condition to the new (face) velocity field. */

  dt = timestep (uf);
}
