extern scalar * tracers;
extern vector u;
extern double dt;

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
