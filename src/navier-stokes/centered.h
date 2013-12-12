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
vector g[];

/**
The parameters are the viscosity coefficient $\mu$ and the specific
volume $\alpha = 1/\rho$ (with default unity). $\alpha$ is a face
vector because we will need its values at the face locations of
the velocity components. 

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

face vector alpha;  // default dt
scalar      alphac; // default dt
face vector mu, a;  // default zero
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
The default velocity and pressure are zero. */

event defaults (i = 0)
{
  const face vector unity[] = {1.,1.};
  alpha = unity;
  const face vector zero[] = {0.,0.};
  a = zero;
  
  CFL = 0.8;
  foreach() {
    foreach_dimension()
      u.x[] = g.x[] = 0.;
    p[] = pf[] = 0.;
  }
  boundary ({p,pf,u,g});
}

/**
We initialise the face velocity field and apply boundary conditions
after user initialisation. */

event init (i = 0)
{
  boundary ({p,u});
  trash ({uf});
  foreach_face()
    uf.x[] = (u.x[] + u.x[-1,0])/2.;
  boundary_normal ({uf});
  boundary_tangent ({uf});
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
  scalar dux[], duy[];
  vector du;
  du.x = dux; du.y = duy;
  if (u.x.gradient) {
    foreach()
      foreach_dimension()
        du.x[] = u.x.gradient (u.x[-1,0], u.x[], u.x[1,0])/Delta;
  }
  else {
    foreach()
      foreach_dimension()
	du.x[] = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
  }
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1,0])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i,0] + g.x[i,0]*dt/2. +
      s*min(1., 1. - s*un)*du.x[i,0]*Delta/2.;
    double fyy = u.y[i,0] < 0. ? u.x[i,1] - u.x[i,0] : u.x[i,0] - u.x[i,-1];
    uf.x[] -= dt*u.y[i,0]*fyy/(2.*Delta);
  }
  boundary ((scalar *){uf});
}

/**
The timestep for this iteration is controlled by the CFL condition
(and the timing of upcoming events). */

event stability (i++,last) {
  dt = dtnext (t, timestep (uf));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. */

event vof (i++,last);
event tracer_advection (i++,last);
event density (i++,last);

event advection_term (i++,last)
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha);
    foreach_dimension()
      advection ({u.x}, uf, dt, g.x); // fixme: list for {g.x}
  }
}

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
  boundary ((scalar *){u});  
}

event viscous_term (i++,last)
{
  if (mu.x) {
    const scalar dtc[] = dt;
    correction (dt);
    mgu = viscosity (u, mu, alphac ? alphac : dtc);
    correction (-dt);
  }

  trash ({uf});
  foreach_face()
    uf.x[] = (u.x[] + u.x[-1,0])/2.;
}

event acceleration (i++,last)
{ 
  if (a.x) {
    (const) face vector af = a;
    boundary_normal ({af});
    foreach_face()
      uf.x[] += dt*af.x[];
    //    foreach_dimension() {
      foreach_boundary (right,false) {
	u.x[] += dt*(af.x[] + af.x[1,0])/2.;
	double ub = (u.x[] + u.x.boundary[right] (point, u.x))/2.;
	u.x[] -= dt*(af.x[] + af.x[1,0])/2. + (uf.x[ghost] - ub)/2.;
	uf.x[ghost] = ub;
      }
      foreach_boundary (left,false) {
	u.x[] += dt*(af.x[] + af.x[1,0])/2.;
	double ub = (u.x[] + u.x.boundary[left] (point, u.x))/2.;
	u.x[] -= dt*(af.x[] + af.x[1,0])/2. + (uf.x[] - ub)/2.;
	uf.x[] = ub;
      }
      foreach_boundary (top,false) {
	u.y[] += dt*(af.y[] + af.y[0,1])/2.;
	double ub = (u.y[] + u.y.boundary[top] (point, u.y))/2.;
	u.y[] -= dt*(af.y[] + af.y[0,1])/2. + (uf.y[ghost] - ub)/2.;
	uf.y[ghost] = ub;
      }
      foreach_boundary (bottom,false) {
	u.y[] += dt*(af.y[] + af.y[0,1])/2.;
	double ub = (u.y[] + u.y.boundary[bottom] (point, u.y))/2.;
	u.y[] -= dt*(af.y[] + af.y[0,1])/2. + (uf.y[] - ub)/2.;
	uf.y[] = ub;
      }
      //    }
  }
  boundary_normal ({uf}); 
}

event projection (i++,last)
{
  mgp = project (uf, p, alpha);
  (const) face vector af = a, alphaf = alpha;
  face vector gf[];
  foreach_face()
    gf.x[] = af.x[] - alphaf.x[]*(p[] - p[-1,0])/(dt*Delta);
  boundary_normal ({gf});
  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1,0])/2.;
  boundary ((scalar *){g});
  // check BCs for g when it is refined
  correction (dt);
}
