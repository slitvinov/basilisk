// Generic predictor/corrector time-integration

#include "utils.h"

// Required from solver
// fields updated by time-integration
extern scalar * evolving;
// how to compute updates
double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;

// User-provided functions
// gradient
double (* gradient)  (double, double, double) = minmod2;

double t = 0., dt = 0.;

trace
static void advance_generic (scalar * output, scalar * input, scalar * updates,
			     double dt)
{
  if (input != output)
    trash (output);
  foreach() {
    scalar o, i, u;
    for (o,i,u in output,input,updates)
      o[] = i[] + dt*u[];
  }
  boundary (output);
}

static void (* advance) (scalar * output, scalar * input, scalar * updates,
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

trace
void run()
{
  t = 0.;
  init_grid (N);

  // main loop
  timer start = timer_start();
  int i = 0; long tnc = 0;
  while (events (i, t, true)) {
    // list of updates
    scalar * updates = list_clone (evolving);
    dt = dtnext (t, update (evolving, updates, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      scalar * predictor = list_clone (evolving);
      /* predictor */
      advance (predictor, evolving, updates, dt/2.);
      /* corrector */
      update (predictor, updates, dt);
      delete (predictor);
      free (predictor);
    }
    advance (evolving, evolving, updates, dt);
    delete (updates);
    free (updates);
    long nc = 0;
    foreach (reduction(+:nc)) 
      nc++;
    tnc += nc;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid();
}
