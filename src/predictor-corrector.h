// Generic predictor/corrector time-integration

#include "utils.h"

// Required from solver
// fields updated by time-integration
extern scalar * evolving;
// how to compute updates
double update (scalar * evolving, scalar * updates, double dtmax);

// User-provided functions
// gradient
double (* gradient)  (double, double, double) = minmod2;

double t = 0., dt = 0.;

static void advance_generic (scalar * output, scalar * input, scalar * updates,
			     double dt)
{
  foreach() {
    scalar o, i, u;
    for (o,i,u in output,input,updates)
      o[] = i[] + dt*u[];
  }
  boundary (output);
}

void (* advance) (scalar * output, scalar * input, scalar * updates,
		  double dt) = advance_generic;

event defaults (i = 0)
{
  // limiting
  for (scalar s in all)
    s.gradient = gradient;

  // default values
  foreach()
    for (scalar s in all)
      s[] = 0.;
  boundary (all);
}

void run()
{
  t = 0.;
  init_grid(N);

  // main loop
  timer start = timer_start();
  int i = 0; long tnc = 0;
  while (events (i, t)) {
    // list of updates
    scalar * updates = clone (evolving);
    dt = dtnext (t, update (evolving, updates, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      /* predictor */
      advance (updates, evolving, updates, dt/2.);
      /* corrector */
      update (updates, updates, dt);
    }
    advance (evolving, evolving, updates, dt);
    delete (updates);
    free (updates);
    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
