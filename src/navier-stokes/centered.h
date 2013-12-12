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

  CFL = 0.8;
  foreach() {
    foreach_dimension()
      u.x[] = 0.;
    p[] = pf[] = 0.;
  }
  boundary ({p,pf,u});
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
  const face vector zero[] = {0.,0.};
  (const) face vector af = a.x ? a : zero;
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1,0])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i,0] + af.x[i,0]*dt/2. +
      s*min(1., 1. - s*un)*g.x[i,0]*Delta/2.;
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

    (const) face vector alphaf = alpha;
    (const) face vector af = a;
    foreach_dimension() {
      scalar dp[];
      foreach()
	dp[] = (af.x[] + af.x[1,0])/2. 
	- (alphaf.x[]*(p[] - p[-1,0]) +
	   alphaf.x[1,0]*(p[1,0] - p[]))/(2.*Delta*dt);
      boundary ({dp});
      advection ({u.x}, uf, dt, dp);
    }
  }
}

static void correction (double dt, face vector a)
{
  double s = sign(dt);
  const face vector zero[] = {0.,0.};
  (const) face vector af = a.x ? a : zero;
  (const) face vector alphaf = alpha;
#if QUADTREE
  face vector g[];
  foreach_face()
    g.x[] = dt*af.x[] - s*alphaf.x[]*(p[] - p[-1,0])/Delta;
  boundary_normal ({g});
  foreach()
    foreach_dimension()
      u.x[] += (g.x[] + g.x[1,0])/2.;
#else
  foreach()
    foreach_dimension()
      u.x[] += dt*(af.x[] + af.x[1,0])/2. -
        s*(alphaf.x[]*(p[] - p[-1,0]) +
	   alphaf.x[1,0]*(p[1,0] - p[]))/(2.*Delta);
#endif
  boundary ((scalar *){u});  
}

event viscous_term (i++,last)
{
  if (mu.x) {
    const scalar dtc[] = dt;
    correction (dt, a);
    mgu = viscosity (u, mu, alphac ? alphac : dtc);
    correction (-dt, a);
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
    foreach()
      foreach_dimension()
        u.x[] += dt*(af.x[] + af.x[1,0])/2.;
    foreach_dimension() {
      foreach_boundary (right,false) {
	double ub = (u.x[] + u.x.boundary[right] (point, u.x))/2.;
	u.x[] -= (uf.x[ghost] - ub)/2.;
	uf.x[ghost] = ub;
      }
      foreach_boundary (left,false) {
	double ub = (u.x[] + u.x.boundary[left] (point, u.x))/2.;
	u.x[] -= (uf.x[] - ub)/2.;
	uf.x[] = ub;
      }
    }
  }
  boundary_normal ({uf}); 
}

event projection (i++,last)
{
  mgp = project (uf, p, alpha);
  correction (1., (vector){0,0});
}
