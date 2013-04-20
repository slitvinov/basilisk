#include "utils.h"

scalar f = new scalar, u = new scalar, v = new scalar;

// user-provided functions
void parameters         (void);
void initial_conditions (void);
void boundary_f         (scalar f);
void boundary_u_v       (scalar u, scalar v);
void boundary_gradient  (scalar fx, scalar fy);

void gradient (const scalar f, vector g)
{
  foreach()
    foreach_dimension()
      g.x[] = (f[1,0] - f[-1,0])/(2.*delta);
  //      g.x[] = minmod ((f[1,0] - f[])/(f[] - f[-1,0]))*(f[] - f[-1,0])/delta;
}

void fluxes_upwind_bcg (const scalar f, const vector g,
			const vector u, 
			vector flux,
			double dt)
{
  foreach()
    foreach_dimension() {
      double un = dt*u.x[]/delta, s = sign(un);
      int i = -(s + 1.)/2.;
      double f2 = f[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*delta/2.;
      double vn = u.y[i,0] + u.y[i,1];
      double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
      f2 -= dt*vn*fyy/(4.*delta);
      flux.x[] = f2*u.x[];
    }
}

double timestep (const scalar u, const scalar v)
{
  double dtmax = DT/CFL;
  foreach(reduction(min:dtmax)) {
    double dx = L0*delta;
    double dt = dx/fabs(u[]);
    if (dt < dtmax) dtmax = dt;
    dt = dx/fabs(v[]);
    if (dt < dtmax) dtmax = dt;
  }
  return dtmax*CFL;
}

void run (void)
{
  parameters ();
  init_grid (N);

  events_init();
  initial_conditions();
  boundary_f(f);

  timer_t start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep (u, v));
    vector flux = new vector, g = new vector, uv = {u,v};
    gradient (f, g);
    boundary_gradient (g.x, g.y);
    fluxes_upwind_bcg (f, g, uv, flux, dt);
    boundary_u_v (flux.x, flux.y);
    foreach()
      f[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/delta;
    boundary_f (f);
    foreach() tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
