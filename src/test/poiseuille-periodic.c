// same as poiseuille.c but with periodic boundary conditions

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

int main() {
  origin (0, -0.5);
  periodic (right);
  
  stokes = true;
  TOLERANCE = 1e-5;

  u.t[top] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
  for (N = 8; N <= 64; N *= 2)
    run(); 
}

scalar un[];

event init (t = 0) {
  // we also check for hydrostatic balance
  const face vector g[] = {1.,1.};
  a = g;
  const face vector muc[] = {1.,1.};
  mu = muc;
  foreach()
    un[] = u.x[];
}

event logfile (t += 0.1; i <= 100) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6)
    return 0; /* stop */
}

event profile (t = end) {
  printf ("\n");
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  scalar e[];
  foreach()
    e[] = u.x[] - 0.5*(0.25 - y*y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
