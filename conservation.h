// generic solver for system of conservation laws
#include "utils.h"

scalar * tendencies = NULL;

// Default user-provided parameters
// conserved quantities
extern scalar * conserved;
// gradient
double (* gradient)    (double, double, double) = minmod2;
// user-provided functions
void      parameters   (void);
void      init         (void);
void      flux         (const double *, double *, double *);

/* generic central-upwind scheme: see e.g. section 3.1 in
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */ 
double riemann (const double * right, const double * left,
		double delta, double * f, int len, 
		double dtmax)
{
  double fm[len], fp[len], em[2], ep[2];
  flux (right, fm, em);
  flux (left,  fp, ep);
  double ap = max(ep[1], em[1]); ap = max(ap, 0.);
  double am = min(ep[0], em[0]); am = min(am, 0.);
  double a = max(ap, -am);
  if (a > 0.) {
    for (int i = 0; i < len; i++)
      f[i] = (ap*fm[i] - am*fp[i] + ap*am*(left[i] - right[i]))/(ap - am);
    double dt = CFL*delta/a;
    if (dt < dtmax)
      dtmax = dt;
  }
  else
    for (int i = 0; i < len; i++)
      f[i] = 0.;
  return dtmax;
}

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
  double r[len], l[len]; // right/left Riemann states
  double f[len];         // fluxes for each conserved quantity

  // compute fluxes and tendencies
  foreach_face() {
    double dx = delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in conserved,slopes) {
      l[i] = s[] - dx*g.x[];
      r[i++] = s[-1,0] + dx*g.x[-1,0];
    }
    // Riemann solver
    dtmax = riemann (r, l, delta, f, len, dtmax);
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
