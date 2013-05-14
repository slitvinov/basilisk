// generic solver for system of conservation laws
#include "utils.h"

typedef struct {
  double l, r;
} state; // left-right Riemann state

scalar * tendencies = NULL;

// Default user-provided parameters
// conserved quantities
extern scalar * conserved;
// gradient
double (* gradient)    (double, double, double) = minmod2;
// user-provided functions
void      parameters   (void);
void      init         (void);
double    riemann      (state * s, double delta, double * flux, double dtmax);

static double fluxes (scalar * conserved, double dtmax)
{
  // allocate space for slopes
  vector * slopes = NULL;
  for (scalar s in conserved) {
    vector slope = new vector;
    foreach_dimension()
      slope.x.gradient = zero; // first-order interpolation only
    slopes = vectors_append (slopes, slope);
  }
  // compute slopes
  gradients (conserved, slopes);

  // allocate space for fluxes
  int len = list_len (conserved);
  state  c[len]; // Riemann states for each conserved quantity
  double f[len]; // fluxes for each conserved quantity

  // compute fluxes and tendencies
  foreach_face() {
    double dx = delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in conserved,slopes) {
      c[i].l = s[] - dx*g.x[];
      c[i++].r = s[-1,0] + dx*g.x[-1,0];
    }
    // Riemann solver
    dtmax = riemann (c, delta, f, dtmax);
    // update tendencies
    i = 0;
    for (scalar ds in tendencies) {
      ds[] += f[i];
      ds[-1,0] -= f[i++];
    }
  }

#if QUADTREE
  // propagate updates from fine to coarse
  foreach_halo()
    for (scalar ds in tendencies) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = 0.;
    }
#endif

  // free space for slopes
  delete ((scalar *) slopes);

  return dtmax;
}

static void update (scalar * conserved1, scalar * conserved2, double dt)
{
  if (conserved1 != conserved2)
    trash (conserved1);
  foreach() {
    scalar s1, s2, ds;
    for (s1,s2,ds in conserved1,conserved2,tendencies) {
      s1[] = s2[] + dt*ds[]/delta;
      ds[] = 0.;
    }
  }
  boundary (conserved1);
}

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

  // allocate tendencies
  tendencies = NULL;
  for (scalar s in conserved) {
    scalar ds = new scalar;
#if QUADTREE
    // we need to reinitialise the tendencies during refinement
    ds.refine = refine_reset;
#endif
    tendencies = list_append (tendencies, ds);
  }

  // limiting
  for (scalar s in conserved)
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
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    dt = dtnext (t, fluxes (conserved, DT));

    if (gradient != zero) {
      /* 2nd-order time-integration */
      // temporary storage
      scalar * conserved1 = clone (conserved);
      /* predictor */
      update (conserved1, conserved, dt/2.);
      /* corrector */
      fluxes (conserved1, dt);
      delete (conserved1);
    }

    update (conserved, conserved, dt);

    foreach (reduction(+:tnc))
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  delete (tendencies);
  free_grid ();
}
