#include "grid/cartesian.h"
#include "advection.h"

void parameters()
{
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
}

double bump (double x, double y)
{
  double r2 = x*x + y*y; 
  double coeff = 20. + 20000.*r2*r2*r2*r2;
  return (1. + cos(20.*x)*cos(20.*y))*exp(-coeff*r2)/2.;
}

void init()
{
  foreach()
    f[] = bump(x,y);
}

#define end 0.785398

// event (t += 0.1; t <= end) output_matrix (f, N, stdout, true);

/* event (t += 0.1; t <= 5.) { */
/*   scalar l = new scalar; */
/*   foreach() l[] = level; */
/*   output_matrix (l, N, stdout, false); */
/* } */

/* event (t += 0.1; t <= 5.) { */
/*   char s[80]; */
/*   FILE * fp; */
/*   scalar l = new scalar; */
/*   foreach() l[] = level; */
/*   sprintf (s, "level-%g", t); */
/*   fp = fopen (s, "w"); */
/*   output_matrix (l, N, fp, false); */
/*   fclose (fp); */

/*   sprintf (s, "f-%g", t); */
/*   fp = fopen (s, "w"); */
/*   output_matrix (f, N, fp, false); */
/*   fclose (fp);   */

/*   sprintf (s, "cells-%g", t); */
/*   fp = fopen (s, "w"); */
/*   output_cells (fp); */
/*   fclose (fp);   */
/* } */

/* event (t += 0.1; t <= 5.) { */
/*   wavelet (f, w); */
/*   output_matrix (w, N, stdout, false); */
/* } */

event velocity (i++) {
  foreach_face(x) u.x[] = -8.*y;
  foreach_face(y) u.y[] =  8.*x;
  boundary_flux (u);
}

event logfile (t = {0,end}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event field (t = end) {
  scalar e[];
  foreach()
    e[] = f[] - bump(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);

  if (N == 256)
    output_field ({e}, stdout, N);
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run ();
}
