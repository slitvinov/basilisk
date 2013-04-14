#include "grid/cartesian1D.h"
#include "saint-venant.h"

void parameters()
{
  //  gradient = generalized_minmod; theta = 1.;
  N = 200;
  CFL = 0.5;
}

void init()
{
  foreach() {
    double epsilon = 1e-5;
    w[] = 1. + epsilon*(-0.4 < x && x < -0.3);
    B[] = 0.25*(cos(pi*x/0.1) + 1.)*(fabs(x) < 0.1);
  }
}

event (t = 0.7) {
  fprintf (stderr, "t: %g dt: %g\n", t, dt);
  foreach()
    printf ("%g %g %.12f %g %g %g\n", x, y, w[], hu[], fhu.x[], fw.x[]);
  printf ("\n");
}

int main() { run(); }
