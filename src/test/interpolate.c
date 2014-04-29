/* interpolation */

#include "utils.h"

scalar v[];

double radius;

int refine_func (Point point, void * data) {
  return x*x + y*y > radius*radius;
}

int main (int argc, char ** argv)
{
  origin (-0.5, -0.5);
  for (int n = 8; n <= 64; n *= 2) {
    init_grid (n);

    radius = 0.49;
    refine_function (refine_func, NULL, NULL);
    radius = 0.55;
    refine_function (refine_func, NULL, NULL);

    foreach()
      v[] = cos(2.*pi*x)*cos(2.*pi*y);
    boundary ({v});

    //    FILE * fp = fopen("error","w");
    double emax = 0.;
    int ni = 4*n + 7;
    double Delta = 1./ni;
    for (int i = 0; i <= ni; i++) {
      double x = Delta*i - 0.5;
      for (int j = 0; j <= ni; j++) {
	double y = Delta*j - 0.5;
	double e = fabs (cos(2.*pi*x)*cos(2.*pi*y) - interpolate (v, x, y));
	//	fprintf (fp, "%g %g %g\n", x, y, e);
	if (e > emax) emax = e;
      }
      //      fprintf (fp, "\n");
    }
    //    fclose (fp);

    fprintf (stderr, "%d %g\n", n*4, emax);
    if (n == 16)
      output_cells (stdout);

    free_grid();
  }
}
