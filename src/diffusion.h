/**
# Time-implicit discretisation of reaction--diffusion equations

We want to discretise implicitly the reaction--diffusion equation
$$ 
\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$ 
where $\beta f + r$ is a reactive term and $D$ is the diffusion
coefficient.

Using a time-implicit backward Euler discretisation, this can be
written
$$
\frac{f^{n+1} - f^{n}}{dt} = \nabla\cdot(D\nabla f^{n+1}) + \beta
f^{n+1} + r
$$
Rearranging the terms we get
$$
\nabla\cdot(D\nabla f^{n+1}) + (\beta - \frac{1}{dt}) f^{n+1} =
- \frac{f^{n}}{dt} - r
$$
This is a Poisson--Helmholtz problem which can be solved with a
multigrid solver. */

#include "poisson.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` is omitted it is set to one. If $\beta$ is omitted it is
set to zero. Both `D` and $\beta$ may be constant fields.

Note that the `r` and $\beta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

struct Diffusion {
  // mandatory
  scalar f, r;
  double dt;
  // optional
  face vector D; // default 1
  scalar beta;        // default 0
};

mgstats diffusion (struct Diffusion p)
{

/**
We define $f$ and $r$ for convenience. The integration time $dt$ is
passed as $- 1/dt$ to fit the parameters of the Poisson--Helmholtz solver. */

  scalar f = p.f, r = p.r;
  double idt = - 1./p.dt;

/**
We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  foreach()
    r[] = idt*f[] - r[];

/**
If $\beta$ is not provided, the diagonal term $\lambda$ is a constant,
otherwise we use $\beta$ to store it. */

  const scalar lambda[] = idt;
  if (p.beta) {
    scalar beta = p.beta;
    foreach()
      beta[] += idt;
    lambda = p.beta;
  }

/**
Finally we solve the system. */

  return poisson (f, r, p.D, lambda);
}
