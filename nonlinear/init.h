#include <stdio.h>
#include <gsl/gsl_integration.h>

/* ---------------- Initial conditions ------------------- */

#define H0 1.
#define FROUDE 0.1

double vtheta (double r) {
  return FROUDE*(r < 0.4)*(1. + cos((r - 0.2)/0.2*M_PI))/2.;
}

double h0p (double r, void * p) {
  double vt = vtheta(r);
  return vt*(2.*OMEGA + vt/r)/G;
}

double h0 (double r) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;
  F.function = &h0p;
  gsl_integration_qags (&F, 0, r, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

void initial_conditions (void * grid)
{
  var h1 = var(h1);
  foreach (grid) {
    val(h1,0,0) = val(h,0,0) = (H0 + h0(sqrt (x*x + y*y)));
    val(u,0,0)  = - vtheta(sqrt (xu*xu + yu*yu))*yu/sqrt (xu*xu + yu*yu);
    val(v,0,0)  =   vtheta(sqrt (xv*xv + yv*yv))*xv/sqrt (xv*xv + yv*yv);
  } end_foreach();
}

/* ------------------ Boundary conditions ------------------- */

void boundary_h (void * grid, var h)
{
  symmetry (grid, h);
}

void boundary_b (void * grid)
{
  symmetry (grid, b);
}

void boundary_ke_psi (void * grid)
{
  symmetry (grid, ke);

  foreach_boundary (grid, top) {
    val(psi,0,1) = (val(v,0,1) - val(v,-1,1) + val(u,0,0) - val(u,0,1))/(L0*delta);
  } end_foreach_boundary();
  foreach_boundary (grid, right) {
    val(psi,1,0) = (val(v,1,0) - val(v,0,0) + val(u,1,-1) - val(u,1,0))/(L0*delta);
  } end_foreach_boundary();
}

void boundary_u (void * grid, var u, var v)
{
  uv_symmetry (grid, u, v);
}

/* ------------------ Output helper functions --------------- */

double error (void * grid)
{
  var e = var(e), h1 = var(h1);
  double max = 0.;
  foreach (grid) {
    val(e,0,0) = fabs (val(h1,0,0)  - val(h,0,0));
    if (val(e,0,0) > max) max = val(e,0,0);
  } end_foreach();
  return max;
}

double energy (void * grid)
{
  double se = 0.;
  foreach (grid)
    se += val(h,0,0)*val(ke,0,0) + G*(val(h,0,0) - H0)*(val(h,0,0) - H0)/2.*delta*delta;
  end_foreach();
  return se*L0*L0;
}

#if 0
void output_field (void * grid, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  int _n = n; /* fixme */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g\n", XC(i), YC(j), e[i][j]);
    fprintf (fp, "\n");
  }
}
#endif
