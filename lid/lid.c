#include "navier-stokes1.h"

void parameters ()
{
  // number of grid points
  N = 64;
  // viscosity
  NU = 1e-3;
  // end time
  TMAX = 300;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

void initial_conditions (void * grid)
{
  /* default to zero */
}

void boundary_p (void * grid, var p, int l)
{
  foreach_boundary_level (grid, right, l)  p[1,0]  = p[];
  foreach_boundary_level (grid, left, l)   p[-1,0] = p[];
  foreach_boundary_level (grid, top, l)    p[0,1]  = p[];
  foreach_boundary_level (grid, bottom, l) p[0,-1] = p[];
}

void boundary_u (void * grid, var u, var v)
{
  foreach_boundary (grid, right) {
    u[1,0] = 0.;
    v[1,0] = - v[];
  }
  foreach_boundary (grid, left) {
    u[-1,0] = - u[1,0];
    u[] = 0.;
    v[-1,0] = - v[];
  }
  foreach_boundary (grid, top) {
    v[0,1] = 0.;
    u[0,1] = 2. - u[];
  }
  foreach_boundary (grid, bottom) {
    v[0,-1] = - v[0,1];
    v[] = 0.;
    u[0,-1] = - u[];
  }
}

static double energy (void * grid)
{
  double se = 0.;
  foreach (grid)
    se += (sq(u[] + u[1,0)] + sq(v[] + v[0,1))]/8.*delta*delta;
  return se*L0*L0;
}

var un = new var; /* we need another variable */

int events (void * grid, int i, double t, double dt)
{
  if (i % 10 == 0) {
    double du = change (grid, u, un);
    if (i > 0 && du < 1e-4)
      return 1; /* stop */
    fprintf (stderr, "t: %f %.9f %g\n", t, energy (grid), du);
  }
  if (i % 100 == 0)
    output_field (grid, u, N, stdout);
  return 0; /* continue */
}

void end (void * grid)
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (grid, u, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (grid, v, x, 0));
  fclose (fp);
}

int main() { run (); }
