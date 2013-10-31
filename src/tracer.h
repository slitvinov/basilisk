/**
# Generic tracer advection event

This event integrates advection equations of the form
$$
\partial_tf_i+\mathbf{u}\cdot\nabla f_i=0
$$
where $\mathbf{u}$ is the velocity field and $f_i$ are a list of
passive tracers.

The `tracers` list is defined elsewhere (typically by the user), the
staggered vector field `u` and the timestep `dt` are defined by a
solver. */

extern scalar * tracers;
extern staggered vector u;
extern double dt;

/**
Given these variables, and a `tracer_fluxes` function, the event
below integrates the advection equation, using an explicit scheme with
timestep `dt`, for each tracer in the list. */

event tracer_advection (i = 1; i++)
{
  for (scalar f in tracers) {
    vector flux[];
    tracer_fluxes (f, u, flux, dt);
    foreach()
      f[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/Delta;
    boundary ({f});
  }
}
