/**
# Tracer advection event

This event integrates advection equations of the form
$$
\partial_tf_i+\mathbf{u_f}\cdot\nabla f_i=0
$$
where $\mathbf{u_f}$ is the velocity field and $f_i$ are a list of
passive tracers.

The `tracers` list is defined elsewhere (typically by the user), the
face vector field `uf` and the timestep `dt` are defined by a
solver. */

extern scalar * tracers;
extern face vector uf;
extern double dt;

/**
The integration is performed using the Bell-Collela-Glaz scheme. */

#include "bcg.h"

event tracer_advection (i++,last) {
  advection (tracers, uf, dt);
}
