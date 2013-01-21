#include <math.h>
#include <stdio.h>

struct _Data {
  double a, da1, da2;
};

#include "cartesian.c"
#include "utils.c"

#define A 0.2

void advance (void * grid, double t, var * f, var * df)
{
  foreach (grid) {
    val(df[0],0,0) = exp(-A*t)*cos(t) - A*val(f[0],0,0);
  } end_foreach();
}

void update (void * grid, double t, var * f) {}

int main ()
{
  void * grid = init_grid(1);
  var a = var(a), da1 = var(da1), da2 = var(da2);
  double dt = 0.4;
  var f[1] = {a}, df[2][1] = {{da1},{da2}};
  int i = 0;
  for (double t = 0.; t <= 6.*3.14159265359; t += dt, i++) {
    foreach (grid) { 
      fprintf (stderr, "%g %g %g\n", t, val(a,0,0), val(a,0,0) - exp(-A*t)*sin(t));
    } end_foreach();
    runge_kutta (grid, t, dt, 2, 1, f, df, advance, update);
  }
  free_grid(grid);
}
