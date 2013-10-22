// generic solver for systems of conservation laws
#include "utils.h"

// user-provided parameters and functions:
// list of conserved scalar quantities
extern scalar * scalars;
// list of conserved vector quantities
extern vector * vectors;
// fluxes/eigenvalues for each conserved quantity
void flux (const double * state, double * flux, double * eigenvalue);

// the list of evolving fields required by predictor-corrector.h
scalar * evolving;

event defaults (i = 0)
{
  // the evolving scalar fields are the conserved scalars + the
  // components of the conserved vectors
  evolving = list_concat (scalars, (scalar *) vectors);
}

event init (i = 0)
{
  boundary (all);
}

/* generic central-upwind scheme: see e.g. section 3.1 in
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */ 
static double riemann (const double * right, const double * left,
		       double Delta, double * f, int len, 
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
    double dt = CFL*Delta/a;
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

  // recover scalars/vectors and the corresponding vector/tensor slopes/fluxes
  int scalars_len = list_len (scalars);

  scalar * scalars = list_copy (conserved);
  if (scalars) scalars[scalars_len] = -1;
  vector * vectors = vectors_from_scalars (&conserved[scalars_len]);

  vector * scalar_slopes = vectors_copy (slopes);
  if (scalar_slopes) scalar_slopes[scalars_len] = (vector){-1,-1};
  tensor * vector_slopes = tensors_from_vectors (&slopes[scalars_len]);

  vector * scalar_fluxes = vectors_copy (lflux);
  if (scalar_fluxes) scalar_fluxes[scalars_len] = (vector){-1,-1};
  tensor * vector_fluxes = tensors_from_vectors (&lflux[scalars_len]);

  // compute fluxes
  foreach_face (reduction (min:dtmax)) {
    // left/right states
    double dx = Delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in scalars,scalar_slopes) {
      l[i] = s[] - dx*g.x[];
      r[i++] = s[-1,0] + dx*g.x[-1,0];
    }
    vector v;
    tensor t;
    for (v,t in vectors,vector_slopes) {
      l[i] = v.x[] - dx*t.x.x[];
      r[i++] = v.x[-1,0] + dx*t.x.x[-1,0];
      l[i] = v.y[] - dx*t.y.x[];
      r[i++] = v.y[-1,0] + dx*t.y.x[-1,0];
    }
    // Riemann solver
    dtmax = riemann (r, l, Delta, f, len, dtmax);
    // update fluxes
    i = 0;
    for (vector fs in scalar_fluxes)
      fs.x[] = f[i++];
    for (tensor fv in vector_fluxes) {
      fv.x.x[] = f[i++];
      fv.y.x[] = f[i++];
    }
  }

  // free space for slopes and lists
  free (scalars);
  free (vectors);
  free (scalar_slopes);
  free (vector_slopes);
  free (scalar_fluxes);
  free (vector_fluxes);
  delete ((scalar *) slopes);
  free (slopes);

  boundary_normal (lflux);

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
      o[] = i[] + dt*(f.x[] - f.x[1,0] + f.y[] - f.y[0,1])/Delta;
  }
  boundary (output);
  delete ((scalar *)lflux);
  lflux = NULL;
}

#include "predictor-corrector.h"
