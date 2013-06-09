// generic solver for system of conservation laws
#include "utils.h"

// user-provided parameters and functions:
// list of conserved quantities
#define evolving conserved
// fluxes/eigenvalues for each conserved quantity
void flux (const double *, double *, double *);
// user-provided initial conditions
void init (void);

void init_internal (void)
{
  init();
  boundary (all);
}

/* generic central-upwind scheme: see e.g. section 3.1 in
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */ 
static double riemann (const double * right, const double * left,
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

// fluxes
vector * lflux = NULL;

double fluxes (scalar * conserved, double dtmax)
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
  for (scalar s in conserved) {
    vector f1 = new vector;
    lflux = vectors_append (lflux, f1);
  }

  // compute fluxes
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
    // update fluxes
    i = 0;
    for (vector f1 in lflux)
      f1.x[] = f[i++];
  }

  // free space for slopes
  delete ((scalar *) slopes);
  free (slopes);

  return dtmax;
}

void update (scalar * output, scalar * input, double dt)
{
  if (output != input)
    trash (output);
  foreach() {
    scalar o, i;
    vector f;
    for (o,i,f in output,input,lflux)
      o[] = i[] + dt*(f.x[] - f.x[1,0])/delta;
  }
  boundary (output);
  delete ((scalar *)lflux);
  lflux = NULL;
}

#include "predictor-corrector.h"
