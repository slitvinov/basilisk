#include <math.h>
#include <stdio.h>
#include "grid/cartesian.h"
#include "utils.h"

#define A 0.2

scalar a = new scalar, da1 = new scalar, da2 = new scalar;

void advance (double t, scalar * f, scalar * df)
{
  foreach()
    val(df[0],0,0) = exp(-A*t)*cos(t) - A*val(f[0],0,0);
}

void update (double t, scalar * f) {}

int main ()
{
  init_grid(1);
  double dt = 0.4;
  scalar f[1] = {a}, df[2][1] = {{da1},{da2}};
  int i = 0;
  for (double t = 0.; t <= 6.*3.14159265359; t += dt, i++) {
    foreach()
      fprintf (stderr, "%g %g %g\n", t, a[], a[] - exp(-A*t)*sin(t));
    runge_kutta (2, t, dt, 1, f, df, advance, update);
  }
  free_grid();
}
