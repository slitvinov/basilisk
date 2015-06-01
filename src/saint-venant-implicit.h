/**
# A semi-implicit Saint-Venant solver

This solver is based on the paper of [Kwatra et al,
2009](references.bib#kwatra2009). The [Saint-Venant
equations](saint-venant.h) are solved semi-implicitly to lift the
timestep restriction due to gravity waves.

We start with the advection solver and add the Poisson solver. */

#include "advection.h"
#include "poisson.h"

/**
The primitive variables are the water depth $h$ and the flux
$\mathbf{q}=h\mathbf{u}$. The acceleration of gravity is *G*. */

vector q[];
scalar h[];
double G = 1.;

/**
Both $\mathbf{q}$ and $h$ are advected using the conservative
advection scheme. */

scalar * tracers = {q,h};

/**
The default slope-limiting is the same as for the explicit
[Saint-Venant solver](saint-venant.h). */

event defaults (i = 0) {
  theta = 1.3;
  gradient = minmod2;
}

/**
The advection solver defines a face velocity field. We initialise it
by interpolation from the flux and the water depth. */

event init (i = 0) {
  boundary ({q,h});
  foreach_face()
    uf.x[] = (q.x[] + q.x[-1])/(h[] + h[-1]);
  boundary ((scalar *){uf});
}

/**
The equation for the pressure is a Poisson--Helmoltz problem which we
will solve with the [multigrid solver](poisson.h). The statistics for
the solver will be stored in *mgp*. */

mgstats mgp;

/**
Up to this point, the advection solver includes events which will
advect $h$ and $\mathbf{q}$ to the next timestep i.e. we now have
$$
\begin{eqnarray}
  h_{n + 1} & = & h_n - \Delta t \nabla \cdot (\mathbf{u}_n h_{n + 1 / 2})\\
  \mathbf{q}_{\star} & = & \mathbf{q}_n - \Delta t \nabla \cdot
  (\mathbf{u}_n \mathbf{q}_{n + 1 / 2})
\end{eqnarray}
$$
We now need to add the pressure gradient term to get
$\mathbf{q}_{n+1}$ from $\mathbf{q}_{\star}$.
*/

event pressure (i++, last)
{

  /**
  We first define a temporary face velocity field $\mathbf{u}_\star$
  using simple averaging from $\mathbf{q}_{\star}$ and $h_{n + 1}$. */

  foreach_face()
    uf.x[] = (q.x[] + q.x[-1])/(h[] + h[-1]);
  boundary ((scalar *){uf});

  /**
  The evolution equation for the pressure is
  $$p_t +\mathbf{u} \cdot \nabla p = - \rho c^2 \nabla \cdot \mathbf{u}$$
  with $\rho$ the density ($h$ in our case) and $c$ the speed of sound
  ($\sqrt{gh}$ for Saint-Venant).

  Following the classical [projection
  method](navier-stokes/centered.h#approximate-projection) for
  incompressible flows, we set
  $$
  \mathbf{u}_{n + 1} = \mathbf{u}_{\star} - 
                         \Delta t \frac{\nabla p}{\rho^{n + 1}}
  $$
  The evolution equation for the pressure can then be discretised as
  $$
  \frac{p_{n + 1} - p_n}{\Delta t} +\mathbf{u}_n \cdot \nabla p_n = 
     - \rho c^2_{n + 1} \nabla \cdot \mathbf{u}_{n + 1}
  $$
  which gives, after some manipulations, the Poisson--Helmholtz equation
  $$
  \lambda_n p_{n + 1} + \nabla \cdot \left( \frac{\nabla p}{\rho_{}}
  \right)_{n + 1} = \lambda_n p_{\star} + \frac{1}{\Delta t} \nabla \cdot
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
  For the specific case of Saint-Venant we have $p=gh^2/2$ and $c^2=gh$. */
  
  scalar p[], lambda[], rhs[];
  foreach() {

    /**
    To compute the right-hand-side of the Poisson--Helmholtz
    equations, we first compute the divergence of the velocity
    field. */
    
    rhs[] = 0.;
    foreach_dimension()
      rhs[] += uf.x[1] - uf.x[];
    rhs[] /= dt*Delta;

    /**
    The pressure $p_\star$ is obtained simply using the equation of
    state ($p=gh^2/2$) and the value of $h$ at time $n+1$. */
    
    p[] = G*sq(h[])/2.;

    /**
    The $\lambda$ coefficient is computed and the pressure term is
    added to the r.h.s. */
    
    lambda[] = -1./(sq(dt)*2.*p[]);
    rhs[] += lambda[]*p[];
  }
  boundary ({p});

  /**
  We then compute the coefficients for the Laplacian operator. */
  
  face vector alpha[];
  foreach_face()
    alpha.x[] = 2./(h[] + h[-1]);
  boundary_flux ({alpha});

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
  We define a face pressure using a density-weighted averaging. */
  
  face vector pf = alpha;
  foreach_face()
    pf.x[] = (p[]*h[] + p[-1]*h[-1])/(h[] + h[-1]);
  boundary_flux ({pf});

  /**
  And finally we apply the pressure gradient term to the flux/momentum. */
  
  foreach()
    foreach_dimension()
      q.x[] -= dt*(pf.x[1] - pf.x[])/Delta;
  boundary ((scalar *){q});
}
