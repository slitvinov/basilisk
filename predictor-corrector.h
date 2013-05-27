// Generic predictor/corrector time-integration

// Required from solver
// fields updated by time-integration
extern scalar * evolving;
// how to compute fluxes
double fluxes (scalar * evolving, double dtmax);
// how to update evolving fields
void   update (scalar * output, scalar * input, double dt);

// User-provided parameters/functions
// gradient
double (* gradient)  (double, double, double) = minmod2;
void      parameters (void);
void      init       (void);

double t = 0., dt = 0.;

void run()
{
  t = 0.;
  parameters();
  init_grid(N);

  // limiting
  for (scalar s in all)
    s.gradient = gradient;

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
  int i = 0; long tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, fluxes (evolving, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      // temporary storage
      scalar * temporary = clone (evolving);
      /* predictor */
      update (temporary, evolving, dt/2.);
      /* corrector */
      fluxes (temporary, dt);
      // free temporary storage
      delete (temporary);
      free (temporary);
    }
    update (evolving, evolving, dt);

    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
