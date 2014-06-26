/**
# A solver for the Green-Naghdi equations

The Green-Naghdi equations (also known as the Serre or fully nonlinear
Boussinesq equations) can be seen as an extension of the [Saint-Venant
equations](saint-venant.h) to the next order in term of the
*shallowness parameter*
$$
\mu = \frac{h_0^2}{L^2}
$$
with $h_0$ the typical depth and $L$ the typical horizontal scale. In
contrast to the Saint-Venant equations the Green-Naghdi equations have
*dispersive* wave solutions.

The solver is built by adding a source term to the momentum equation
of the [Saint-Venant solver](saint-venant.h). Following [Bonneton et
al, 2011](/src/references.bib#bonneton2011), this source term can be
written in one dimension
$$
\partial_t \left( hu \right) + \ldots = h \left( \frac{g}{\alpha}
   \partial_x \eta - D \right)
$$
where $D$ verifies
$$
\alpha h\mathcal{T} \left( D \right) + hD = b
$$
and
$$
b = h \left( \frac{g}{\alpha} \partial_x \eta + 2 h \partial_x h \left(
\partial_x u \right)^2 + \frac{4}{3} h^2 \partial_x u \partial_x^2 u \right)
$$
This linear system can be inverted with the multigrid Poisson
solver. We declare *D* as a global variable so that it can be re-used
as initial guess for the Poisson solution. The solver statistics will
be stored in *mgD*. */

#include "saint-venant.h"
#include "poisson.h"

scalar D[];
mgstats mgD;

/**
To add the source term to the Saint-Venant system we overload the
default *sources* function with this one. The function takes
a list of the current evolving scalar fields and fills the
corresponding *updates* with the source terms. */

static void green_naghdi (scalar * current, scalar * updates)
{
  scalar h = current[0];
  vector u = { current[1], current[2] };
  double alpha = 1.;

  /**
  We first compute the right-hand-side $b$. */

  scalar b[];
  foreach() {
    double dxu = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    b[] = h[]*(G/alpha*(eta[1,0] - eta[-1,0])/2. +
	       h[]*(h[1,0] - h[-1,0])*sq(dxu) +
	       4./3.*sq(h[])*dxu*(u.x[1,0] + u.x[-1,0] - 2.*u.x[])/Delta)/Delta;
  }

  /**
  The equation for $D$ can be rewritten (using the definition of
  $\mathcal{T}$)
  $$
  - \partial_x \left( \alpha \frac{h^3}{3} \partial_x D \right) + hD = 
  \partial_x \left( \beta \partial_x D \right) + hD = b
  $$
  This is a Poisson-Helmholtz problem which we solve with the
  multigrid solver. */

  face vector beta[];
  foreach_face() {
    double hf = (h[] + h[-1,0])/2.;
    beta.x[] = - alpha*hf*hf*hf/3.;
  }
  boundary_normal ({beta});
  mgD = poisson (D, b, beta, h);

  /**
  We can then compute the updates for $h$ (0) and $hu$. Note that we
  need to be careful here as *current* and *updates* can be identical
  i.e. *h* and *dh*, *u* and *dhu* can be identical. */

  scalar dh = updates[0];
  vector dhu = { updates[1], updates[2] };

  foreach() {
    dhu.x[] = h[]*(G/alpha*(eta[1,0] - eta[-1,0])/(2.*Delta) - D[]);
    dh[] = dhu.y[] = 0.;
  }
}

event defaults (i = 0) {
  sources = green_naghdi;
  foreach()
    D[] = 0.;
  boundary ({D});
}
