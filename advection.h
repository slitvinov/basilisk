#include "utils.h"

scalar f[]; // tracer
vector u[]; // velocity

// Default parameters
// gradient
double (* gradient) (double, double, double) = NULL; // centered
// user-provided functions
void parameters  (void);
void init        (void);

void fluxes_upwind_bcg (const scalar f, const vector g,
			const vector u, 
			vector flux,
			double dt)
{
  foreach_face() {
    double un = dt*u.x[]/delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*delta/2.;
    double vn = u.y[i,0] + u.y[i,1];
    double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
    f2 -= dt*vn*fyy/(4.*delta);
    flux.x[] = f2*u.x[];
  }
  boundary_flux (flux);
}

double timestep (const vector u)
{
  double dtmax = DT/CFL;
  foreach(reduction(min:dtmax)) {
    double dx = L0*delta;
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = dx/fabs(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  }
  return dtmax*CFL;
}

void run (void)
{
  parameters();
  init_grid (N);
  f.gradient = gradient;
  foreach()
    f[] = u.x[] = u.y[] = 0.;
  init();
  boundary ({f, u.x, u.y});

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep (u));
    vector flux[], g[];
    gradients ({f}, {g});
    fluxes_upwind_bcg (f, g, u, flux, dt);
    foreach()
      f[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/delta;
    boundary ({f});
    foreach (reduction (+:tnc))
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
