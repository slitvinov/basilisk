// similar to stream.c but using the centered NS solver
#include "navier-stokes/centered.h"

#define MAXLEVEL 8

void parameters()
{
  X0 = Y0 = -0.5;
  N = 1 << MAXLEVEL;
}

event init (t = 0)
{
  scalar psi[], omega[];

  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

  double dd = 0.1;
  foreach() {
    omega[] = (exp(-(sq(x - dd) + sq(y))/(dd/10.)) +
	       exp(-(sq(x + dd) + sq(y))/(dd/10.)));
    psi[] = 0.;
  }
  boundary ({psi,omega});
  poisson (psi, omega);
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});
}

void vorticity (vector u, scalar omega)
{
  foreach()
    omega[] = (u.y[1,0] - u.y[-1,0] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
  boundary ({omega});
}

event logfile (t = {0,30}) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);
  fprintf (stderr, "%g %d %g %g %d\n", t, i, dt, s.sum, mgp.i);
}

event movie (t += 0.2; t <= 30) {
  static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, fp, linear = true);

  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  foreach()
    omega[] = level;
  output_ppm (omega, fp1);
}

event output (t += 5) {
  static int nf = 0;
  scalar omega[];
  vorticity (u, omega);
  printf ("file: omega-%d\n", nf);
  output_field ({omega}, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({l});
  nf++;
}

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-4,1e-4},
		 MAXLEVEL, list = {p,pf,u,uf});
}
#endif

int main() { run(); }
