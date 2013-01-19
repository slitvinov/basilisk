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

double ** h1, ** e;

void initial_conditions (void * m, int n)
{
  h1 = matrix_new (n, n, sizeof (double));
  e = matrix_new (n, n, sizeof (double));
  foreach (m, n) {
    h1[I][J] = h(0,0) = (H0 + h0(sqrt (x*x + y*y)));
    u(0,0) = - vtheta(sqrt (xu*xu + yu*yu))*yu/sqrt (xu*xu + yu*yu);
    v(0,0) = vtheta(sqrt (xv*xv + yv*yv))*xv/sqrt (xv*xv + yv*yv);
  } end_foreach();
}

/* ------------------ Boundary conditions ------------------- */

void boundary_h (void * m, int n, var h)
{
  symmetry (m, n, h);
}

void boundary_b (void * m, int n)
{
  symmetry (m, n, var(b));
}

void boundary_ke_psi (void * m, int n)
{
  symmetry (m, n, var(ke));

  double dx = L0/n;
  foreach_boundary (m, n, top) {
    psi(0,1) = (v(0,1) - v(-1,1) + u(0,0) - u(0,1))/dx;    
  } end_foreach_boundary();
  foreach_boundary (m, n, right) {
    psi(1,0) = (v(1,0) - v(0,0) + u(1,-1) - u(1,0))/dx;
  } end_foreach_boundary();
}

void boundary_u (void * m, int n, var u, var v)
{
  uv_symmetry (m, n, u, v);
}

/* ------------------ Output helper functions --------------- */

double error (void * m, int n)
{
  double max = 0.;
  foreach (m, n) {
    e[I][J] = fabs (h1[I][J]  - h(0,0));
    if (e[I][J] > max) max = e[I][J];
  } end_foreach();
  return max;
}

double energy (void * m, int n)
{
  double se = 0.;
  foreach (m, n)
    se += h(0,0)*ke(0,0) + G*(h(0,0) - H0)*(h(0,0) - H0)/2.;
  end_foreach();
  return se*(L0/n)*(L0/n); /* fixme */
}

void output_field (void * m, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  int _n = n; /* fixme */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g\n", XC(i), YC(j), e[i][j]);
    fprintf (fp, "\n");
  }
}
