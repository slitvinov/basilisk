#define X0 -0.5
#define Y0 -0.5
#define XC (L0*(i - 0.5)/n + X0)
#define YC (L0*(j - 0.5)/n + Y0)
#define XU (L0*(i - 1.)/n + X0)
#define YU (L0*(j - 0.5)/n + Y0)
#define XV (L0*(i - 0.5)/n + X0)
#define YV (L0*(j - 1.)/n + Y0)

double H0 = 1., ETA0 = 0.1, R0 = 0.1;

void initial_conditions (double ** u, double ** v, double ** h, double ** b, int n)
{
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      u[i][j] = v[i][j] = 0.;
      double x = XC, y = YC;
#if 0
      b[i][j] = 0.8*exp(-5.*(x - 0.9)*(x - 0.9) - 50.*(y - 0.5)*(y - 0.5));
      h[i][j] = (1. + 0.01*exp(-2000.*(x - 0.1)*(x - 0.1)) - b[i][j]);
#else
      b[i][j] = 0.;
      h[i][j] = (H0 + ETA0*exp (-(x*x + y*y)/(R0*R0)));
      x = XU; y = YU;
      u[i][j] = 2.*G*ETA0*y/(F0*R0*R0)*exp (-(x*x + y*y)/(R0*R0));
      x = XV; y = YV;
      v[i][j] = - 2.*G*ETA0*x/(F0*R0*R0)*exp (-(x*x + y*y)/(R0*R0));
#endif
    }
  symmetry_conditions (h, n);
  solid_walls_conditions (u, v, n);
}

static double error (double ** h, int n)
{
  double max = 0.;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      double x = XC, y = YC;
      double e = fabs ((H0 + ETA0*exp (-(x*x + y*y)/(R0*R0))) - h[i][j]);
      if (e > max) max = e;
    }
  return max;
}
