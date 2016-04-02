/**
# Generic time loop

The `run()` function below implements a generic time loop which
executes events until termination. */

#include "utils.h"

/**
The time `t` and timestep `dt` can be accessed as global variables. */

double t = 0., dt = 1.;

trace
void run (void)
{
  t = 0.; dt = 1.;
  init_grid (N);

  int i = 0;
  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (i, t, true)) {

    /**
    We store the total number of cells advanced in time for computing
    speed statistics. */

    update_perf();
    i++; t = tnext;
  }

  /**
  Time/speed statistics are written out on standard output. */

  timer_print (perf.gt, i, perf.tnc);

  free_grid();
}
