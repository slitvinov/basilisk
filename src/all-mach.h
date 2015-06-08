/**
# An "all Mach" flow solver

We wish to solve the generic momentum equation
$$
\partial_t\mathbf{q} + \nabla\cdot(\mathbf{q}\mathbf{u}) = 
- \nabla p + \mathbf{a}
$$
with $\mathbf{q}=\rho\mathbf{u}$ the momentum, $\mathbf{u}$ the
velocity, $\rho$ the density, $p$ the pressure and $\mathbf{a}$ an
acceleration. The pressure is defined through an equation of state and
verifies the evolution equation
$$
\partial_t p + \mathbf{u}\cdot\nabla p = -\rho c^2\nabla\cdot\mathbf{u}
$$
with $c$ the speed of sound. 

We build the solver using the generic time loop, Bell-Collela-Glaz
advection scheme and multigrid Poisson--Helmholtz solver. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "poisson.h"

/**
The primitive variables are the momentum $\mathbf{q}$, pressure $p$
and (face) velocity field $\mathbf{uf}$. */

vector q[];
scalar p[];
face vector uf[];

/**
The equation of state is defined by the pressure field *ps* and $\rho
c^2$. Both fields are zero by default (i.e. the fluid is
incompressible). The reciprocal volume $\alpha=1/\rho$ is one by
default (i.e. the density is one). */

scalar ps[];
(const) scalar rhoc2 = zeroc;
(const) face vector alpha = unityf, a = zerof;

/**
We need a field to store the combined pressure gradient and
acceleration field. */

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
  
#if QUADTREE
  for (scalar s in {q}) {
    s.refine = s.prolongation = refine_linear;
    s.coarsen = coarsen_volume_average;
  }
#endif
}

event init (i = 0) {
  boundary ({q,p});

  /**
  The face velocity field is obtained by simple linear interpolation
  from the momentum field. We make sure that $\alpha$ is defined by
  calling the "properties" event (see below). */
  
  event ("properties");
  foreach_face()
    uf.x[] = alpha.x[]*(q.x[] + q.x[-1])/2.;
  boundary ((scalar *){uf});
}

/**
The timestep is computed using the CFL condition on the face velocity
field. */

event stability (i++,last) {
  dt = dtnext (t, timestep (uf));
}

/**
Tracers are advected by these events. */

event vof (i++,last);
event tracer_advection (i++,last);

/**
Momentum is advected using the [BCG scheme](bcg.h) and including the
combined acceleration and pressure gradient source term *g* i.e. we now have
$$
\mathbf{q}_\star = \mathbf{q}_n - \Delta t
\nabla\cdot(\mathbf{q}_{n+1/2}\mathbf{u}) 
$$
*/

event advection_term (i++,last) {
  advection ((scalar *){q}, uf, dt, (scalar *){g});
}

/**
The equation of state (i.e. fields $\alpha$, $\rho c^2$ and *ps*) is
defined by this event. */

event properties (i++,last)
{
  boundary_flux ({alpha});

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
  using simple averaging from $\mathbf{q}_{\star}$ and the
  acceleration term. */

  foreach_face()
    uf.x[] = alpha.x[]*(q.x[] + q.x[-1])/2. + dt*a.x[];
  boundary ((scalar *){uf});
  
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
  face velocity field at time $n + 1$. */
  
  foreach_face()
    uf.x[] -= dt*alpha.x[]*(p[] - p[-1])/Delta;
  boundary ((scalar *){uf});
  
  /**
  We define a combined face pressure and acceleration using a
  density-weighted averaging. */
  
  face vector gf[];
  foreach_face()
    gf.x[] = a.x[] - alpha.x[]*(p[] - p[-1])/Delta;
  boundary_flux ({gf});

  /**
  We use the face pressure pressure to get the centered pressure gradient. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[]/alpha.x[] + gf.x[1]/alpha.x[1])/2.;
  boundary ((scalar *){g});

  /**
  And finally we apply the pressure gradient term to the flux/momentum. */
  
  foreach()
    foreach_dimension()
      q.x[] += dt*g.x[];
  boundary ((scalar *){q});
}
