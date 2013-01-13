#define X0 -0.5
#define Y0 -0.5
#define XC (L0*(i - 0.5)/n + X0)
#define YC (L0*(j - 0.5)/n + Y0)
#define XU (L0*(i - 1.)/n + X0)
#define YU (L0*(j - 0.5)/n + Y0)
#define XV (L0*(i - 0.5)/n + X0)
#define YV (L0*(j - 1.)/n + Y0)

#include <gsl/gsl_integration.h>

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

void initial_conditions (double ** u, double ** v, double ** h, double ** b, int n)
{
  h1 = matrix_new (n + 2, n + 2, sizeof (double));
  e = matrix_new (n + 2, n + 2, sizeof (double));
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      double x = XC, y = YC;
      b[i][j] = 0.;
      h1[i][j] = h[i][j] = (H0 + h0(sqrt (x*x + y*y)));
      x = XU; y = YU;
      u[i][j] = - vtheta(sqrt (x*x + y*y))*y/sqrt (x*x + y*y);
      x = XV; y = YV;
      v[i][j] = vtheta(sqrt (x*x + y*y))*x/sqrt (x*x + y*y);
    }
  symmetry_conditions (h, n);
  solid_walls_conditions (u, v, n);
}

static double error (double ** h, int n)
{
  double max = 0., maxh = 0.;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      e[i][j] = fabs (h1[i][j]  - h[i][j]);
      if (e[i][j] > max) max = e[i][j];
      if (fabs (h1[i][j]) > maxh) maxh = fabs (h1[i][j]);
    }
  return max/maxh;
}

static double energy (double ** u, double ** v, double ** h, int n)
{
  double se = 0.;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      se += h[i][j]*KE(u,v,i,j) + G*(h[i][j] - H0)*(h[i][j] - H0)/2.;
  return se*(L0/n)*(L0/n);
}
