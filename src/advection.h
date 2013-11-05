/**
# An advection solver

We wish to solve the advection equations
$$
\partial_tf_i+\mathbf{u}\cdot\nabla f_i=0
$$
where $\mathbf{u}$ is the velocity field and $f_i$ are a list of
passive tracers.  This can be done with a flux-based advection scheme
such as the 2nd-order, unsplit, upwind scheme of [Bell-Collela-Glaz,
1989](references.bib#bell89).

The main time loop is defined in [run.h](). A stable timestep needs to
respect the [CFL
condition](http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition). The
Bell-Collela-Glaz flux computation is defined in [bcg.h](). */

#include "run.h"
#include "timestep.h"
#include "bcg.h"

/**
We allocate the (staggered) velocity field. The `gradient` function is
used to set the type of slope-limiting required. The default is to not
use any limiting (i.e. a purely centered slope estimation). */

staggered vector u[]; // velocity

// Default parameters
double (* gradient) (double, double, double) = NULL; // centered

/**
The default boundary conditions are symmetry on all boundaries
i.e. the (staggered) normal components of velocity are zero on all
boundaries. */

u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

/**
Here we set the default staggering for velocity and the gradient
functions for each tracer (as defined in the user-provided `tracers`
list). We also set default values (zero) for the tracer and velocity
fields. */

event defaults (i = 0)
{
  for (scalar f in tracers)
    f.gradient = gradient;
  foreach()
    for (scalar f in tracers)
      f[] = 0.;
  boundary (tracers);
  foreach_face()
    u.x[] = 0.;
  boundary ((scalar *){u});
}

/**
The timestep is set using the velocity field and the CFL
criterion. This is done after all other events (the `last`
keyword). The integration itself is performed in the events of
[tracer.h]() (included from [bcg.h]()). */

event velocity (i++, last) {
  dt = timestep (u);
}
