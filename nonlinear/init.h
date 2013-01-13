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

void initial_conditions (Data * m, int n)
{
  h1 = matrix_new (n + 2, n + 2, sizeof (double));
  e = matrix_new (n + 2, n + 2, sizeof (double));
  foreach (m) {
    double x = XC, y = YC;
    b(0,0) = 0.;
    h1[i][j] = hn(0,0) = (H0 + h0(sqrt (x*x + y*y)));
    x = XU; y = YU;
    un(0,0) = - vtheta(sqrt (x*x + y*y))*y/sqrt (x*x + y*y);
    x = XV; y = YV;
    vn(0,0) = vtheta(sqrt (x*x + y*y))*x/sqrt (x*x + y*y);
  }
}

static double error (Data * m, int n)
{
  double max = 0.;
  foreach (m) {
    e[i][j] = fabs (h1[i][j]  - h(0,0));
    if (e[i][j] > max) max = e[i][j];
  }
  return max;
}

static double energy (Data * m, int n)
{
  double se = 0.;
  foreach (m)
    se += h(0,0)*KE(m,i,j,n) + G*(h(0,0) - H0)*(h(0,0) - H0)/2.;
  return se*(L0/n)*(L0/n);
}
