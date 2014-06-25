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

This solver is built following the ideas of [Le Métayer et al,
2010](/src/references.bib#lemetayer2010).

Using the non-local change of variable
$$
K = \left( 1 +\mathcal{T} \right) u
$$
with the operator $\mathcal{T}$ defined as
$$
\mathcal{T}w = - \frac{1}{3 h} \nabla \left( h^3 \nabla\cdot w \right)
$$
the Green-Naghdi system of equations can be written (in one dimension)
$$
\partial_th + \partial_x(hu) = 0
$$
$$
\partial_thK + \partial_x\left(hKu +\frac{gh^2}{2}+\alpha\right) = 0
$$
with
$$
\alpha = -\frac{2}{3}h^3(\partial_xu)^2
$$
This system can be solved using the solver for systems of conservation
laws, and the operator $\mathcal{T}$ can be inverted with the Poisson
solver (to recover $u$ from $K$). */

#include "conservation.h"
#include "poisson.h"

/**
The primary variables are the fluid depth $h$ and the velocity $u$. We
need to declare $hK$ and $\alpha$ as conserved variables so that they
are accessible when computing fluxes. The Poisson solver statistics
are stored in *mgu*. */

scalar u[], h[], hK[], alpha[];
scalar * scalars = {h, hK, u, alpha};
vector * vectors = NULL;

mgstats mgu;
double G = 1.;

/**
The fluxes for the conservative variables $h$ and $hK$ are computed as
functions of the reconstructed values of $h$, $hK$, $u$ and
$\alpha$. The fluxes for $u$ and $\alpha$ are set to zero. */

void flux (const double * s, double * f, double e[2])
{
  double h = s[0], hK = s[1], u = s[2], alpha = s[3];
  f[0] = h*u;
  f[1] = hK*u + G*h*h/2. + alpha;
  f[2] = f[3] = 0.;
  // min/max eigenvalues
  double c = sqrt(G*h);
  e[0] = u - c; // min
  e[1] = u + c; // max
}

/**
## Diagnostic equations for $u$ and $\alpha$

Although we declared $u$ and $\alpha$ as conserved variables, they are
not, rather we obtain their updated values through diagnostic
equations linking them with the conserved variables $h$ and $hK$. */

static void u_from_hK (scalar u, scalar hK, scalar h, scalar alpha)
{
  
  /**
  Knowing $h$ and $hK$, we can compute the value of $u$ by inverting
  the operator $\mathcal{T}$. In one dimension, we rearrange the terms
  to obtain
  $$
  \frac{\partial}{\partial x}
  \left(- \frac{h^3}{3}\frac{\partial u}{\partial x}\right) + hu = hK
  $$
  This is a Poisson-Helmholtz equation which can be solved with the
  [generic Poisson
  solver](poisson.h#application-to-the-poissonhelmholtz-equation). */

  face vector beta[];
  foreach_face() {
    double hf = (h[] + h[-1,0])/2.;
    beta.x[] = - hf*hf*hf/3.;
  }
  boundary_normal ({beta});
  mgu = poisson (u, hK, beta, h);

  /**
  Once we have $u$, we can compute $\alpha$ using simple centered
  finite differences. */

  foreach() {
    double dudx = (u[1,0] - u[-1,0])/(2.*Delta);
    alpha[] = - 2./3.*h[]*h[]*h[]*sq(dudx);
  }
  boundary ({alpha});
}

/**
## Time update

The default update function of the predictor-corrector scheme computes
the new values by simply adding the computed conservative fluxes to
the current values. We also need to update the values of $u$ and
$\alpha$ (using the function above). To do this, we overload the
default update function. */

static void advance_green_naghdi (scalar * output, scalar * input, 
				  scalar * updates, double dt)
{
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  scalar hKi = input[1], hKo = output[1], dhK = updates[1];
  scalar ui = input[2], uo = output[2], alphao = output[3];

  // new fields in ho[], hKo[]
  foreach() {
    ho[] = hi[] + dt*dh[];
    hKo[] = hKi[] + dt*dhK[];
    uo[] = ui[]; // we use ui as initial guess for the Poisson solver
  }
  boundary ({ho, hKo, uo});
  
  u_from_hK (uo, hKo, ho, alphao);
}

event defaults (i = 0)
  advance = advance_green_naghdi;

/**
## Initial conditions

The user should only work with $h$ and $u$. We use the definition of
$K$ to compute initial conditions for $hK$ from the initial conditions
for $h$ and $u$ (given by the user). Note that this is just the direct
application of operator $1 + \mathcal{T}$ to $u$ i.e. the reciprocal
function of *u_from_hK()* above.

Finally, we use *u_from_hK()* to obtain consistent initial conditions
for $u$ and $\alpha$. */

event init (i = 0) {
  face vector h3u[];
  foreach_face() {
    double hf = (h[] + h[-1,0])/2.;
    h3u.x[] = hf*hf*hf*(u[] - u[-1,0])/Delta;
  }
  boundary_normal ({h3u});
  foreach() {
    hK[] = h[]*u[] - (h3u.x[1,0] - h3u.x[])/(3.*Delta);
    u[] = 0.;
  }
  boundary ({hK,u});

  u_from_hK (u, hK, h, alpha);
}
