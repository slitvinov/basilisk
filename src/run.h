/**
# Generic time loop

The `run()` function below implements a generic time loop which
executes events until termination. */

#include "utils.h"

/**
The `parameters()` function is called before grid initialisation. It
is usually provided by the user. 

The time `t` and timestep `dt` can be accessed as global variables. */

void parameters (void);
double t = 0., dt = 0.;

void run (void)
{
  t = dt = 0.;
  parameters();
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
