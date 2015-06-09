/**
# An "all Mach" flow solver

We wish to solve the generic momentum equation
$$
\partial_t\mathbf{q} + \nabla\cdot(\mathbf{q}\mathbf{u}) = 
- \nabla p + \rho\mathbf{a}
$$
with $\mathbf{q}=\rho\mathbf{u}$ the momentum, $\mathbf{u}$ the
velocity, $\rho$ the density, $p$ the pressure and $\mathbf{a}$ an
acceleration. The pressure is defined through an equation of state and
verifies the evolution equation
$$
\partial_t p + \mathbf{u}\cdot\nabla p = -\rho c^2\nabla\cdot\mathbf{u}
$$
with $c$ the speed of sound. By default the solver sets $c=\infty$,
$\rho=1$ and the pressure equation reduces to 
$$
\nabla\cdot\mathbf{u} = 0
$$

The advection of momentum is not performed by this solver (so that
different schemes can be used) i.e. in the end, by default, we solve
the incompressible (linearised) Euler equations with a projection
method.

We build the solver using the generic time loop and the multigrid
Poisson--Helmholtz solver. */

#include "run.h"
#include "timestep.h"
#include "poisson.h"

/**
The primitive variables are the momentum $\mathbf{q}$, pressure $p$,
density $\rho$ and (face) velocity field $\mathbf{uf}$. */

vector q[];
scalar p[];
face vector uf[];
(const) scalar rho = unity;

/**
The equation of state is defined by the pressure field *ps* and $\rho
c^2$. Both fields are zero by default (i.e. the fluid is
incompressible). */

scalar ps[];
(const) scalar rhoc2 = zeroc;
(const) face vector a = zerof;

/**
We store the combined pressure gradient and acceleration field in
*g*. */

vector g[];

/**
We initialise default values for the primitive fields. */

event defaults (i = 0) {
  foreach() {
    p[] = 0.;
    foreach_dimension()
      q.x[] = g.x[] = 0.;
  }
  boundary ({q,p,g});
}

event init (i = 0) {
  boundary ({q,p,rho});

  /**
  The face velocity field is obtained by simple linear interpolation
  from the momentum field. We make sure that the density field is
  defined by calling the "properties" event (see below). */
  
  event ("properties");
  foreach_face()
    uf.x[] = (q.x[] + q.x[-1])/(rho[] + rho[-1]);
  boundary ((scalar *){uf});
}

/**
The timestep is computed using the CFL condition on the face velocity
field. */

event stability (i++,last) {
  dt = dtnext (t, timestep (uf));
}

/**
Tracers (including momentum $\mathbf{u}$) are advected by these events. */

event vof (i++,last);
event tracer_advection (i++,last);

/**
The equation of state (i.e. fields $\rho$, $\rho c^2$ and *ps*) is
defined by this event. */

event properties (i++,last)
{

  /**
  If the acceleration is not constant, we reset it to zero. */
  
  if (!is_constant(a.x)) {
    face vector af = a;
    foreach_face()
      af.x[] = 0.;
  }
}

/**
This event can be overloaded to add acceleration terms. */

event acceleration (i++, last)
{
  boundary_flux ({a});
}

/**
The equation for the pressure is a Poisson--Helmoltz problem which we
will solve with the [multigrid solver](poisson.h). The statistics for
the solver will be stored in *mgp*. */

mgstats mgp;

event pressure (i++, last)
{

  /**
  We first define a temporary face velocity field $\mathbf{u}_\star$
  using simple averaging from $\mathbf{q}_{\star}$, $\rho_{n+1}$ and
  the acceleration term. */

  foreach_face()
    uf.x[] = (q.x[] + q.x[-1])/(rho[] + rho[-1]) + dt*a.x[];
  boundary ((scalar *){uf});

  /**
  We compute the specific volume $\alpha=1/\rho$ which will be used in
  the Poisson--Helmholtz equation for the pressure. */
  
  face vector alpha, gf[];
  if (is_constant(rho)) {
    const face vector alphac[] = {1./constant(rho),
				  1./constant(rho),
				  1./constant(rho)};
    alpha = alphac;
  }
  else {
    alpha = gf;
    foreach_face()
      alpha.x[] = 2./(rho[] + rho[-1]);
    boundary_flux ({alpha});
  }
    
  /**
  The evolution equation for the pressure is
  $$p_t +\mathbf{u} \cdot \nabla p = - \rho c^2 \nabla \cdot \mathbf{u}$$
  with $\rho$ the density and $c$ the speed of sound. Following the
  classical [projection
  method](navier-stokes/centered.h#approximate-projection) for
  incompressible flows, we set
  $$
  \mathbf{u}_{n + 1} = \mathbf{u}_{\star} - \Delta t (\alpha\nabla p)_{n+1}
  $$
  The evolution equation for the pressure can then be discretised as
  $$
  \frac{p_{n + 1} - p_n}{\Delta t} +\mathbf{u}_n \cdot \nabla p_n = 
     - \rho c^2_{n + 1} \nabla \cdot \mathbf{u}_{n + 1}
  $$
  which gives, after some manipulations, the Poisson--Helmholtz equation
  $$
  \lambda_{n + 1} p_{n + 1} + \nabla \cdot \left( \alpha \nabla p
  \right)_{n + 1} = \lambda_{n + 1} p_{\star} + \frac{1}{\Delta t} \nabla \cdot
  \mathbf{u}_{\star}
  $$
  with
  $$
  p_{\star} = p_n - \Delta t\mathbf{u}_n \cdot \nabla p_n
  $$
  and
  $$
  \lambda = \frac{- 1}{\Delta t^2 \rho c^2}
  $$
  */

  scalar lambda = rhoc2, rhs = ps;
  foreach() {

    /**
    We compute $\lambda$ and the corresponding term in the
    right-hand-side of the Poisson--Helmholtz equation. */

    if (constant(lambda) == 0.)
      rhs[] = 0.;
    else {
      lambda[] = -1./(sq(dt)*rhoc2[]);
      rhs[] = lambda[]*ps[];
    }
      
    /**
    We add the divergence of the velocity field to the right-hand-side. */

    double div = 0.;
    foreach_dimension()
      div += uf.x[1] - uf.x[];
    rhs[] += div/(dt*Delta);
  }
  
  /**
  The Poisson--Helmholtz solver is called with a [definition of the
  tolerance](poisson.h#377) identical to that used for incompressible
  flows. */
  
  mgp = poisson (p, rhs, alpha, lambda, tolerance = TOLERANCE/sq(dt));

  /**
  The pressure gradient is applied to $\mathbf{u}_\star$ to obtain the
  face velocity field at time $n + 1$. 
  
  We also store the combined face pressure gradient and acceleration
  field *gf*. */

  foreach_face() {
    double dp = 2./(rho[] + rho[-1])*(p[] - p[-1])/Delta;
    uf.x[] -= dt*dp;
    gf.x[] = a.x[] - dp;
  }
  boundary_flux ({gf});

  /**
  And finally we apply the pressure gradient/acceleration term to the
  flux/momentum. We also store the centered, combined pressure
  gradient and acceleration field *g*. */
  
  foreach()
    foreach_dimension() {
      g.x[] = rho[]*(gf.x[] + gf.x[1])/2.;
      q.x[] += dt*g.x[];
    }
  boundary ((scalar *){q,g,uf});
}
