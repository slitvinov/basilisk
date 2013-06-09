#include "grid/multigrid.h"
#include "navier-stokes1.h"

void parameters()
{ 
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // number of grid points
  N = 64;
  // viscosity
  NU = 1e-3;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

u.x[top]    = 2. - u.x[]; // u.x[] = 1. on top boundary
// no slip walls on all other boundaries */
u.x[bottom] = - u.x[];
u.y[left]   = - u.y[];
u.y[right]  = - u.y[];

void init() {}

static double energy()
{
  double se = 0.;
  foreach(reduction(+:se))
    se += (sq(u.x[] + u.x[1,0)] + sq(u.y[] + u.y[0,1))]/8.*delta*delta;
  return se*L0*L0;
}

scalar un[]; /* we need another scalar */

event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event logfile (i += 10; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-4)
    return 1; /* stop */
  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

event outputfile (i += 100) output_matrix (u.x, stdout, N, linear = true);

void end()
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}

int main() { run (); }
