#include <math.h>
#include <stdio.h>
#include "grid/cartesian.h"
#include "utils.h"

#define A 0.2

var a = new var, da1 = new var, da2 = new var;

void advance (void * grid, double t, var * f, var * df)
{
  foreach (grid)
    val(df[0],0,0) = exp(-A*t)*cos(t) - A*val(f[0],0,0);
}

void update (void * grid, double t, var * f) {}

int main ()
{
  void * grid = init_grid(1);
  double dt = 0.4;
  var f[1] = {a}, df[2][1] = {{da1},{da2}};
  int i = 0;
  for (double t = 0.; t <= 6.*3.14159265359; t += dt, i++) {
    foreach (grid)
      fprintf (stderr, "%g %g %g\n", t, a(0,0), a(0,0) - exp(-A*t)*sin(t));
    runge_kutta (2, grid, t, dt, 1, f, df, advance, update);
  }
  free_grid(grid);
}
