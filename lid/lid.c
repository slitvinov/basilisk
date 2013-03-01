#include "navier-stokes1.h"

void parameters()
{
  // number of grid points
  N = 64;
  // viscosity
  NU = 1e-3;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

void initial_conditions()
{
  /* default to zero */
}

void boundary_p (scalar p, int l)
{
  foreach_boundary_level (right, l)  p[1,0]  = p[];
  foreach_boundary_level (left, l)   p[-1,0] = p[];
  foreach_boundary_level (top, l)    p[0,1]  = p[];
  foreach_boundary_level (bottom, l) p[0,-1] = p[];
}

void boundary_u (scalar u, scalar v)
{
  foreach_boundary (right) {
    u[1,0] = 0.;
    v[1,0] = - v[];
  }
  foreach_boundary (left) {
    u[-1,0] = - u[1,0];
    u[] = 0.;
    v[-1,0] = - v[];
  }
  foreach_boundary (top) {
    v[0,1] = 0.;
    u[0,1] = 2. - u[];
  }
  foreach_boundary (bottom) {
    v[0,-1] = - v[0,1];
    v[] = 0.;
    u[0,-1] = - u[];
  }
}

static double energy ()
{
  double se = 0.;
  foreach(reduction(+:se))
    se += (sq(u[] + u[1,0)] + sq(v[] + v[0,1))]/8.*delta*delta;
  return se*L0*L0;
}

scalar un = new scalar; /* we need another scalar */

event (i += 10) {
  double du = change (u, un);
  if (i > 0 && du < 1e-4)
    return 1; /* stop */
  fprintf (stderr, "t: %f %.9f %g\n", t, energy(), du);
}

event (i += 100) output_field (u, N, stdout);

void end()
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (v, x, 0));
  fclose (fp);
}

int main() { run (); }
