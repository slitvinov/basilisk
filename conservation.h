// generic solver for system of conservation laws
#include "utils.h"

// user-provided parameters and functions:
// list of conserved quantities
#define evolving conserved
// fluxes/eigenvalues for each conserved quantity
void flux (const double *, double *, double *);

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

double tendencies (scalar * evolution, scalar * conserved, double dtmax)
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
    for (scalar ds in evolution) {
      ds[] += f[i];
      ds[-1,0] -= f[i++];
    }
  }

#if QUADTREE
  // propagate tendencies from fine to coarse
  foreach_halo()
    for (scalar ds in evolution) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = 0.;
    }
#endif

  // free space for slopes
  delete ((scalar *) slopes);
  free (slopes);

  return dtmax;
}

void update (scalar * conserved1, 
	     scalar * conserved2, scalar * evolution, double dt)
{
  if (conserved1 != conserved2)
    trash (conserved1);
  foreach() {
    scalar s1, s2, ds;
    for (s1,s2,ds in conserved1,conserved2,evolution) {
      s1[] = s2[] + dt*ds[]/delta;
      ds[] = 0.;
    }
  }
  boundary (conserved1);
}

#include "predictor-corrector.h"

