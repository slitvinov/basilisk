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

void init()
{
  u[top] = 2. - u[]; // u[] = 1. on top boundary
  v[left] = - v[];   // no slip walls on all other boundaries
  v[right] = - v[];
  u[bottom] = - u[];
}

static double energy()
{
  double se = 0.;
  foreach(reduction(+:se))
    se += (sq(u[] + u[1,0)] + sq(v[] + v[0,1))]/8.*delta*delta;
  return se*L0*L0;
}

scalar un = new scalar; /* we need another scalar */

int event (i += 10) {
  double du = change (u, un);
  if (i > 0 && du < 1e-4)
    return 1; /* stop */
  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

int event (i += 100) output_matrix (u, N, stdout, true);

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
