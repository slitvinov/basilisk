#include "utils.h"

// Default parameters
void parameters  (void); // user-provided
// time and timestep
double t = 0., dt = 0.;

void run (void)
{
  t = dt = 0.;
  parameters();
  init_grid (N);

  timer start = timer_start();
  int i = 0; long tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, dt);
    foreach(reduction(+:tnc))
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid();
}
