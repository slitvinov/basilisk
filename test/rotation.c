#if ADAPT
# include "grid/quadtree.h"
# include "wavelet.h"
# include "adapt.h"
#else
# include "grid/cartesian.h"
#endif
#include "advection.h"

void boundary_f (scalar f)
{
  boundary (f);
#if ADAPT
  restriction (f, f);
  update_halo (-1, f, f);
#endif
}

void boundary_u_v (scalar u, scalar v)
{
  boundary_uv (u, v);
#if ADAPT
  restriction_u_v (u, v);
  update_halo_u_v (-1, u, v);
#endif
}

void boundary_gradient (scalar fx, scalar fy)
{
#if ADAPT
  restriction (fx, fy);
  update_halo (-1, fx, fy);
#endif
}

void parameters ()
{
  // CFL number
  CFL = 0.8;
}

double bump (double x, double y)
{
  double r2 = x*x + y*y; 
  double coeff = 20. + 20000.*r2*r2*r2*r2;
  return (1. + cos(20.*x)*cos(20.*y))*exp(-coeff*r2)/2.;
}

void initial_conditions ()
{
  foreach()
    f[] = bump(x,y);
}

event (i++) {
#if ADAPT
  scalar w = new scalar;
  restriction (f, f);
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
/*   restriction (f, f); */
/*   wavelet (f, w); */
/*   output_matrix (w, N, stdout, false); */
/* } */

event (i++) {
  foreach() {
    u[] = -8.*yu;
    v[] =  8.*xv;
  }
  boundary_u_v (u, v);
}

event (t = {0,end}) {
  double sum = 0., min = 1e100, max = -1e100;
  foreach(reduction(+:sum) reduction(max:max) reduction(min:min)) {
    sum += f[]*delta*delta;
    if (f[] > max) max = f[];
    if (f[] < min) min = f[];
  }
  fprintf (stderr, "# %f %.12f %g %g\n", t, sum, min, max);  
}

event (t = end) {
  double max = 0., norm1 = 0., norm2 = 0., area = 0.;
  scalar e = new scalar;
  foreach(reduction(max:max) reduction(+:norm1) reduction(+:norm2) reduction(+:area)) {
    e[] = f[] - bump(x,y);
    if (fabs(e[]) > max) max = fabs(e[]);
    double a = sq(delta);
    norm1 += a*fabs(e[]);
    norm2 += a*e[]*e[];
    area  += a;
  }
  fprintf (stderr, "%d %g %g %g\n", N, norm1/area, sqrt(norm2/area), max);
  
  if (N == 256)
    output_matrix (e, N, stdout, false);
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run ();
}
