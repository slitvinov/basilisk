/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(\mu\nabla\mathbf{u})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2012](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "viscosity.h"

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

face vector alpha;  // default one
scalar      alphac; // default one
face vector mu, a;  // default zero
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
## Initial conditions

The default acceleration is zero. The default density is unity.  The
default velocity and pressure are zero. */

event defaults (i = 0)
{
  const face vector zero[] = {0.,0.};
  a = zero;
  const face vector unity[] = {1.,1.};
  alpha = unity;
  
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

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event stability (i++,last) {
  dt = dtnext (t, timestep (uf));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);

/**
The density field (at time $t+\Delta t/2$) can be defined by
overloading this event. */

event density (i++,last);

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz advection scheme and the pressure gradient
and acceleration terms at time $t$ (stored in vector $\mathbf{g}$). */

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
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++,last)
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt);
    advection ((scalar *){u}, uf, dt, (scalar *){g});
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
  boundary ((scalar *){u});  
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will replaced by their values at time
$t+\Delta t$. */

event viscous_term (i++,last)
{
  if (mu.x) {
    const scalar unity[] = 1.;
    correction (dt);
    mgu = viscosity (u, mu, alphac ? alphac : unity, dt);
    correction (-dt);
  }

/**
The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by simple interpolation. */

  trash ({uf});
  foreach_face()
    uf.x[] = (u.x[] + u.x[-1,0])/2.;
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. */

event acceleration (i++,last)
{

/**
We first reset the centered acceleration/pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = 0.;

/**
We then add the acceleration term to the face velocity field. */

  if (a.x) {
    (const) face vector af = a;
    boundary_normal ({af});
    foreach_face()
      uf.x[] += dt*af.x[];

/**
We then apply the boundary conditions to the provisionary face
velocity field.

When Dirichlet conditions are applied to the normal velocity, the
Neumann condition applied on the pressure field is not necessarily
consistent with the cancellation of the acceleration term. While this
is not an issue for the face velocity field (whose value is specified
on the boundary), this is a problem for the centered velocity field:
as the pressure gradient applied to the centered velocity field is the
average of the face pressure gradients, the normal component of the
centered velocity field close to a boundary will only see half the
pressure gradient required to balance the acceleration term. To
compensate for this, we substract from the centered acceleration term
$\mathbf{g}$ half the acceleration term, but only in the case where
Dirichlet conditions are applied to the normal component. */

    foreach_dimension() {
      foreach_boundary (right,false) {
	double ac = dt*(af.x[] + af.x[1,0])/2.;
	u.x[] += ac;
	double ub = (u.x[] + u.x.boundary[right] (point, u.x))/2.;
	u.x[] -= ac;
	g.x[] -= (uf.x[ghost] - ub)/(2.*dt);
	uf.x[ghost] = ub;
      }
      foreach_boundary (left,false) {
	double ac = dt*(af.x[] + af.x[1,0])/2.;
	u.x[] += ac;
	double ub = (u.x[] + u.x.boundary[left] (point, u.x))/2.;
	u.x[] -= ac;
	g.x[] -= (uf.x[] - ub)/(2.*dt);
	uf.x[] = ub;
      }
    }
  }
  boundary_normal ({uf}); 
}

/**
## Approximate projection

To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). */

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt);

/**
We then compute a face field $\mathbf{g}_f$ combining both
acceleration and pressure gradient. */

  (const) face vector af = a, alphaf = alpha;
  face vector gf[];
  foreach_face()
    gf.x[] = af.x[] - alphaf.x[]*(p[] - p[-1,0])/Delta;
  boundary_normal ({gf});

/**
We average these face values to obtain the centered, combined
acceleration and pressure gradient field. */

  foreach()
    foreach_dimension()
      g.x[] += (gf.x[] + gf.x[1,0])/2.;
  boundary ((scalar *){g});

/**
And finally add this term to the centered velocity field. */

  // check BCs for g when it is refined
  correction (dt);
}