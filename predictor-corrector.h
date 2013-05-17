// Generic predictor/corrector time-integration

// Required from solver
// fields updated by time-integration
extern scalar * evolving;
// how to compute tendencies
double tendencies (scalar * evolution, scalar * evolving, double dtmax);
// how to update evolving fields from tendencies
void update (scalar * evolving1, 
	     scalar * evolving, scalar * evolution, double dt);

// User-provided parameters/functions
// gradient
double (* gradient)  (double, double, double) = minmod2;
void      parameters (void);
void      init       (void);

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

  // limiting
  for (scalar s in all)
    s.gradient = gradient;

  // allocate tendencies
  scalar * evolution = NULL;
  for (scalar s in evolving) {
    scalar t = new scalar;
    evolution = list_append (evolution, t);
  }
#if QUADTREE
  // we need the tendencies to be reinitialised during refinement
  for (scalar ds in evolution)
    ds.refine = refine_reset;
#endif

  // default values
  foreach()
    for (scalar s in all)
      s[] = 0.;
  boundary (all);

  // user-defined initial conditions
  init();
  boundary (all);

  // main loop
  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, tendencies (evolution, evolving, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      // temporary storage
      scalar * temporary = clone (evolving);
      /* predictor */
      update (temporary, evolving, evolution, dt/2.);
      /* corrector */
      tendencies (evolution, temporary, dt);
      // free temporary storage
      delete (temporary);
      free (temporary);
    }
    update (evolving, evolving, evolution, dt);

    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free (evolution);
  free_grid ();
}
