#include "grid/cartesian.h"
#include "advection.h"

int refine_right (Point point, void * data)
{
  return x < -0.1;
}

int refine_circle (Point point, void * data)
{
  return x*x + y*y < 0.25*0.25;
}

void parameters()
{
  // coordinates of lower-left corner
  X0 = Y0 = -0.5;
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
}

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

void init()
{
  //  refine_function (0, -1, refine_circle, NULL);
  //  flag_halo_cells ();

  foreach()
    f[] = bump(x,y);
}

int event (i++) {
#if ADAPT
  scalar w = new scalar;
  wavelet (f, w);

  double cmax = 1e-3;
  int nf = refine_wavelet (f, f, w, cmax);
  int nc = coarsen_wavelet (w, cmax/4.);
  flag_halo_cells ();
  boundary_f (f);
#if 0
  fprintf (stderr, "%d refined %d cells coarsened %d cells\n", i, nf, nc);
  check_two_one();
#endif
#endif
}

/* event (t += 0.1; t <= 5.) output_matrix (f, N, stdout, true); */

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
/*   restriction (f, f); */
/*   wavelet (f, w); */
/*   output_matrix (w, N, stdout, false); */
/* } */

int event (i++) {
  foreach_face(x)
    u.x[] = 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*cos((y + 0.5)*pi);
  foreach_face(y)
    u.y[] = - 1.5*sin(2.*pi*t/5.)*cos((x + 0.5)*pi)*sin((y + 0.5)*pi);
  boundary_flux (u);
}

int event (t = {0,5}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

int event (t = 5) {
  scalar e = new scalar;
  foreach()
    e[] = f[] - bump(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  
  if (N == 256)
    output_field (scalars (e), N, stdout, false);
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run ();
}
