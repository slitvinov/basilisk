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
2011](/src/references.bib#lagree2011)).

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

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

p[right] = neumann(a.n[ghost]*fm.n[ghost]/alpha.n[ghost]);
p[left]  = neumann(-a.n[]*fm.n[]/alpha.n[]);

#if AXI
uf.n[bottom] = 0.;
#else // !AXI
#  if dimension > 1
p[top]    = neumann(a.n[ghost]*fm.n[ghost]/alpha.n[ghost]);
p[bottom] = neumann(-a.n[]*fm.n[]/alpha.n[]);
#  endif
#  if dimension > 2
p[front]  = neumann(a.n[ghost]*fm.n[ghost]/alpha.n[ghost]);
p[back]   = neumann(-a.n[]*fm.n[]/alpha.n[]);
#  endif
#endif

/**
## Initial conditions */

event defaults (i = 0)
{

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
    boundary ((scalar *){alpha});
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;
#endif
}

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
  boundary ((scalar *){uf});

  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last) {
  boundary ({alpha, mu, rho});
}

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension()
        du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
  else
    foreach()
      foreach_dimension()
        du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    #endif
    #if dimension > 2
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    #endif
    uf.x[] *= fm.x[];
  }
  boundary ((scalar *){uf});

  delete ((scalar *){du});
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
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
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
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    correction (-dt);
  }

  /**
  The (provisionary) face velocity field at time $t+\Delta t$ is
  obtained by simple interpolation. We also reset the acceleration
  field (if it is not a constant). */

  face vector af = a;
  trash ({uf,af});
  foreach_face() {
    uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
    if (!is_constant(af.x))
      af.x[] = 0.;
  }
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
  boundary ((scalar *){a});
  foreach_face()
    uf.x[] += dt*fm.x[]*a.x[];
  boundary ((scalar *){uf});
}

/**
## Approximate projection

To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). */

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);

  /**
  We then compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = a.x[] - alpha.x[]/fm.x[]*(p[] - p[-1])/Delta;
  boundary_flux ({gf});

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/2.;
  boundary ((scalar *){g});

  /**
  And finally add this term to the centered velocity field. */

  correction (dt);
}

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. */

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif
