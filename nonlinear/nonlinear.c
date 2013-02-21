#include <gsl/gsl_integration.h>
#include "atmosphere.h"

void parameters() {
  // number of grid points
  N = 128;
  // size of the box
  L0 = 1.;
  // acceleration of gravity
  G = 1.;
  // Coriolis parameter
  F0 = 2.*OMEGA;
  // Viscosity
  //  NU = 7e-5;
  NU = 0.;
  // end time
  TMAX = 5.;
  // CFL number
  //  CFL = 0.5;
  CFL = 0.25;
}

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

var h1 = new var;

void initial_conditions (void * grid)
{
  foreach (grid) {
    h1[] = h[] = (H0 + h0(sqrt (x*x + y*y)));
    u[] = - vtheta(sqrt (xu*xu + yu*yu))*yu/sqrt (xu*xu + yu*yu);
    v[] =   vtheta(sqrt (xv*xv + yv*yv))*xv/sqrt (xv*xv + yv*yv);
  }
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

void boundary_u (void * grid, var u, var v)
{
  uv_symmetry (grid, u, v);
}

/* ------------------ Output helper functions --------------- */

var e = new var;

double error (void * grid)
{
  double max = 0.;
  foreach (grid) {
    e[] = fabs (h1[]  - h[]);
    if (e[] > max) max = e[];
  }
  return max;
}

double energy (void * grid)
{
  double se = 0.;
  foreach (grid)
    se += (h[]*ke[] + G*h[]*(h[]/2. + b[]))*delta*delta;
  return se*L0*L0;
}

int events (void * grid, int i, double t, double dt)
{
  if (i % 10 == 0)
    fprintf (stderr, "t: %g %g %g\n", t, error (grid), energy (grid));
  //  if (i % 100 == 0)
  //    output_field (grid, e, stdout);
  return 0;
}

int main() { run(); }