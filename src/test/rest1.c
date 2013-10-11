// lake at rest with emerged island
#include "saint-venant.h"

void parameters()
{
  X0 = Y0 = -0.5;
  N = 32;
}

int refine (Point point, void * data)
{
  return x < 0.1 && y < 0.1;
}

event init (i = 0)
{
  refine_function (refine, NULL, all);

  foreach() {
    zb[] = 0.3*exp(-100*(x*x + y*y));
    h[] = max (0.1 - zb[], 0.);
  }
}

event logfile (i = 1)
{
  norm n = normf (u.x);
  fprintf (stderr, "# %.10f %.10f %.10f\n", n.avg, n.rms, n.max);
  printf ("x y h zb u.x u.y eta\n");
  foreach ()
    printf ("%g %g %g %g %.3g %.3g %.3g\n", x, y, h[], zb[], 
	    u.x[] < 1e-10 ? 0. : u.x[], 
	    u.y[] < 1e-10 ? 0. : u.y[], 
	    h[] > dry ? h[] + zb[] - 0.1 : undefined);
}

int main() { run(); }
