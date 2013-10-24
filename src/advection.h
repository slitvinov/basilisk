#include "utils.h"
#include "timestep.h"
#include "bcg.h"

vector u[]; // velocity

// Default parameters
double (* gradient) (double, double, double) = NULL; // centered
void parameters  (void); // user-provided

// default is no flow through boundaries
u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

// timestep
double dt = 0.;

event defaults (i = 0)
{
  u.x.d.x = u.y.d.y = -1; // staggering for u.x, u.y
  for (scalar f in tracers)
    f.gradient = gradient;
  foreach()
    for (scalar f in tracers)
      f[] = 0.;
  boundary (tracers);
  foreach_face()
    u.x[] = 0.;
  boundary_mac ({u});
}

void run (void)
{
  parameters();
  init_grid (N);

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, timestep (u));
    foreach(reduction(+:tnc))
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
