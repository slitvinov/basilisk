#include <stdio.h>
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

void * matrix_new (int n, int p, int size)
{
  int i;
  void ** m;
  char * a;
  
  m = malloc (n*sizeof (void **));
  a = malloc (n*p*size);
  for (i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return (void *) m;
}

void matrix_free (void * m)
{
  free (((void **) m)[0]);
  free (m);
}

void initial_conditions (Data * m, int n)
{
  h1 = matrix_new (n, n, sizeof (double));
  e = matrix_new (n, n, sizeof (double));
  foreach (m, n) {
    double x = XC(I), y = YC(J);
    b(0,0) = 0.;
    h1[I][J] = hn(0,0) = (H0 + h0(sqrt (x*x + y*y)));
    x = XU(I); y = YU(J);
    un(0,0) = - vtheta(sqrt (x*x + y*y))*y/sqrt (x*x + y*y);
    x = XV(I); y = YV(J);
    vn(0,0) = vtheta(sqrt (x*x + y*y))*x/sqrt (x*x + y*y);
  } end_foreach();
}

static double error (Data * m, int n)
{
  double max = 0.;
  foreach (m, n) {
    e[I][J] = fabs (h1[I][J]  - h(0,0));
    if (e[I][J] > max) max = e[I][J];
  } end_foreach();
  return max;
}

static double energy (Data * m, int n)
{
  double se = 0.;
  foreach (m, n)
    se += h(0,0)*ke(0,0) + G*(h(0,0) - H0)*(h(0,0) - H0)/2.;
  end_foreach();
  return se*(L0/n)*(L0/n); /* fixme */
}

void output_field (Data * m, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  int _n = n; /* fixme */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g\n", XC(i), YC(j), e[i][j]);
    fprintf (fp, "\n");
  }
}
