/**
# Generic time loop

The `run()` function below implements a generic time loop which
executes events until termination. */

#include "utils.h"

/**
The time `t` and timestep `dt` can be accessed as global variables. */

double t = 0., dt = 1.;

void run (void)
{
  t = 0.; dt = 1.;
  init_grid (N);

  timer start = timer_start();
  int i = 0; long tnc = 0;
  while (events (i, t)) {

    /**
    We store the total number of cells advanced in time for computing
    speed statistics. */

    foreach(reduction(+:tnc))
      tnc++;
    i++; t = tnext;
  }

  /**
  Time/speed statistics are written out on standard output. */

  timer_print (start, i, tnc);

  free_grid();
}
