/* interpolation */

#include "utils.h"

scalar v = new scalar;

int main (int argc, char ** argv)
{
  for (int n = 32; n <= 256; n *= 2) {
    init_grid (n);

    foreach()
      v[] = cos(2.*pi*x)*cos(2.*pi*y);
    boundary (v);

    double emax = 0.;
    int ni = n + 7;
    double delta = 1./ni;
    for (int i = 0; i <= ni; i++) {
      double x = delta*i - 0.5;
      for (int j = 0; j <= ni; j++) {
	double y = delta*j - 0.5;
	double e = fabs (cos(2.*pi*x)*cos(2.*pi*y) - interpolate (v, x, y));
	if (e > emax) emax = e;
      }
    }
    fprintf (stderr, "%d %g\n", n, emax);
    free_grid();
  }
}
