#include "utils.h"
#include "bcg.h"
#include "timestep.h"

scalar f[]; // tracer
vector u[]; // velocity

// Default parameters
// gradient
double (* gradient) (double, double, double) = NULL; // centered
// user-provided functions
void parameters  (void);
void init        (void);

// default is no flow through boundaries
u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

// timestep
double dt = 0.;

void run (void)
{
  parameters();
  init_grid (N);

  u.x.d.x = u.y.d.y = -1; // staggering for u.x, u.y
  f.gradient = gradient;
  foreach_face()
    u.x[] = 0.;
  foreach()
    f[] = 0.;
  init();
  boundary ({f});
  boundary_mac ({u});

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, timestep (u));
    vector flux[];
    fluxes_upwind_bcg (f, u, flux, dt);
    foreach(reduction(+:tnc)) {
      f[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/Delta;
      tnc++;
    }
    boundary ({f});
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
