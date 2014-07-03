/**
# A solver for the Green-Naghdi equations

The Green-Naghdi equations (also known as the Serre or fully nonlinear
Boussinesq equations) can be seen as an extension of the [Saint-Venant
equations](saint-venant.h) to the next order $O(\mu^2)$ in term of the
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
written
$$
\partial_t \left( hu \right) + \ldots = h \left( \frac{g}{\alpha}
   \nabla \eta - D \right)
$$
where $D$ verifies
$$
\alpha h\mathcal{T} \left( D \right) + hD = b
$$
and
$$
b = \left[ \frac{g}{\alpha} \nabla \eta +\mathcal{Q}_1 \left( u \right)
\right]
$$
With $\mathcal{T}$ a linear operator to be defined below, as well as
$\mathcal{Q}_1 \left( u \right)$.

This linear system can be inverted with the multigrid Poisson
solver. We declare *D* as a global variable so that it can be re-used
as initial guess for the Poisson solution. The solver statistics will
be stored in *mgD*. */

#include "saint-venant.h"
#include "poisson.h"

scalar D[];
mgstats mgD;

/**
We first define some useful macros, following the notations in
[Bonneton et al, 2011](/src/references.bib#bonneton2011).

Simple centered-difference approximations of the first and second
derivatives of a field. */

#define dx(s)  ((s[1,0] - s[-1,0])/(2.*Delta))
#define d2x(s) ((s[1,0] + s[-1,0] - 2.*s[])/sq(Delta))

/**
The definitions of the $\mathcal{R}_1$ and $\mathcal{R}_2$ operators
$$
\begin{array}{lll}
  \mathcal{R}_1 \left[ h, z_b \right] w & = & - \frac{1}{3 h} \nabla \left(
  h^3 w \right) - \frac{h}{2} w \nabla z_b\\
  & = & - h \left[ \frac{h^{}}{3} \nabla w + w \left( \nabla h + \frac{1}{2}
  \nabla z_b \right)\right]\\
  \mathcal{R}_2 \left[ h, z_b \right] w & = & \frac{1}{2 h} \nabla \left( h^2
  w \right) + w \nabla z_b\\
  & = & \frac{h}{2} \nabla w + w \nabla \left( z_b + h \right)
\end{array}
$$ */

#define R1(h,zb,w) (-h[]*(h[]/3.*dx(w) + w[]*(dx(h) + dx(zb)/2.)))
#define R2(h,zb,w) (h[]/2.*dx(w) + w[]*(dx(zb) + dx(h)))

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
  We first compute the right-hand-side $b$. The general form for the
  $\mathcal{Q}_1$ operator is (eq. (9) of Bonneton et al, 2011).
  $$
  \mathcal{Q}_1 \left[ h, z_b \right] \left( V \right) = -
  2\mathcal{R}_1 \left( \partial_1 V \cdot \partial_2 V^{\perp} + \left(
  \nabla \cdot V \right)^2 \right) +\mathcal{R}_2 \left( V \cdot \left(
  V \cdot \nabla \right) \nabla z_b \right)
  $$
  In one-dimension, this gives
  $$
  \mathcal{Q}_1 \left[ h, z_b \right](u) = - 2\mathcal{R}_1(c) +\mathcal{R}_2(d)
  $$
  with $c = (\partial_xu_x)^2$ and $d = u_x^2\partial^2_xz_b$.
  */

  scalar b[], c[], d[], lambda[];
  foreach() {
    c[] = sq(dx(u.x));
    d[] = sq(u.x[])*d2x(zb);
  }
  boundary ({c,d});

  /**
  The general form for $\mathcal{T}$ is
  $$
  \mathcal{T} \left[ h, z_b \right] w = \mathcal{R}_1 \left[ h, z_b
  \right] \left( \nabla \cdot w \right) +\mathcal{R}_2 \left[ h, z_b \right]
  \left( \nabla z_b \cdot w \right)
  $$
  which gives in one dimension the Poisson-Helmholtz problem
  $$
  - \partial_x \left( \alpha \frac{h^3}{3} \partial_x D \right) + h
  \left[ \alpha \left( \partial_x \eta \partial_x z_b + \frac{h}{2}
  \partial^2_x z_b \right) + 1\right] D = \partial_x \left( \beta
  \partial_x D \right) + \lambda D = b
  $$
  We can now compute $b$ and $\lambda$. We also make sure that
  $\lambda$ is non-zero by setting a minimal value for $h$ (*dry* is a
  small number defined in the Saint-Venant solver) so that the system
  is invertible even in dry areas. */

  foreach() {
    double dxzb = dx(zb), dxeta = dxzb + dx(h);
    b[] = h[]*(G/alpha*dxeta - 2.*R1(h,zb,c) + R2(h,zb,d));
    lambda[] = max(h[],dry)*(alpha*(dxeta*dxzb + h[]/2.*d2x(zb)) + 1.);
  }

  /**
  Finally we solve the Poisson-Helmholtz problem with the multigrid
  solver. */

  face vector beta[];
  foreach_face() {
    double hf = (h[] + h[-1,0])/2.;
    beta.x[] = - alpha*hf*hf*hf/3.;
  }
  boundary_normal ({beta});
  mgD = poisson (D, b, beta, lambda);

  /**
  We can then compute the updates for $h$ (zero) and $hu$. Note that we
  need to be careful here as *current* and *updates* can be identical
  i.e. *h* and *dh*, *u* and *dhu* can be identical. */

  scalar dh = updates[0];
  vector dhu = { updates[1], updates[2] };

  foreach() {

    /**
    We only apply the Green-Naghdi source term when the slope of the
    free surface is smaller than one. The idea is to turn off
    dispersion in areas where the wave is "breaking" (i.e. hydraulic
    jumps). */

    if (fabs(dx(eta)) < 1.)
      dhu.x[] = h[]*(G/alpha*dx(eta) - D[]);
    else
      dhu.x[] = 0.;
    dh[] = dhu.y[] = 0.;
  }
}

/**
In the default setup, we replace the default source terms with our
function. */

event defaults (i = 0) {
  sources = green_naghdi;
  foreach()
    D[] = 0.;
  boundary ({D});
}
