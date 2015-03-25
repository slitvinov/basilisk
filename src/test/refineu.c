/* tangential interpolation on face vector fields  */

#include "utils.h"

scalar h[];
face vector u[];

static double solution (double x, double y)
{
#if 0
  double R0 = 0.1;
  return exp(-(x*x+y*y)/sq(R0));
#else
  return x - y;
#endif
}

// fixme: this does not work with -DTRASH=1

u.n[right]  = solution(x,y);
u.n[left]   = solution(x,y);
u.n[top]    = solution(x,y);
u.n[bottom] = solution(x,y);

u.t[top]    = dirichlet(solution(x,y));
u.t[bottom] = dirichlet(solution(x,y));
u.t[right]  = dirichlet(solution(x,y));
u.t[left]   = dirichlet(solution(x,y));

double tolerance = 1e-4;

static int error()
{
  double max = 0, maxv = 0.;

  scalar eu[];
  foreach(reduction(max:max) reduction(max:maxv))
    for (int i = 0; i <= 1; i++) {
      double xu = x + (i - 0.5)*Delta, yu = y;
      eu[] = fabs (solution(xu,yu) - u.x[i,0]);
      if (eu[] > max)
	max = eu[];

      double xv = x, yv = y + (i - 0.5)*Delta;
      double e = solution(xv,yv) - u.y[0,i];
      if (fabs(e) > maxv)
	maxv = fabs(e);
#if DEBUG
      printf ("%g %g %d %d %g %g\n", 
	      xu, yu, level, cell.neighbors, u.x[i,0], eu[]);
#endif
    }

  fprintf (ferr, "maximum error: %g %g\n", max, maxv);
  stats s = statsf (eu);
  fprintf (ferr, "eu: avg: %g stddev: %g max: %g\n", 
	   s.sum/s.area, s.stddev, s.max);

  return (max != maxv);
}

int main (int argc, char ** argv)
{
  int n = 1024;
  init_grid (n);

  origin (-0.5, -0.5);
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/sq(R0));
  boundary ({h});

  foreach_face()
    u.x[] = solution(x,y);
  boundary ((scalar *){u});

  astats s = adapt_wavelet ({h}, &tolerance, 10);
  while (s.nc) {
    fprintf (ferr, "refined: %d coarsened: %d\n", s.nf, s.nc);
    s = adapt_wavelet ({h}, &tolerance, 10);
  }
  
#if DEBUG
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);
#endif

  if (error())
    return 1;

  scalar div[];
  foreach()
    div[] = u.x[] - u.x[1,0] + u.y[] - u.y[0,1];
  stats sdiv = statsf(div);
  fprintf (ferr, "div before: %g\n", sdiv.max);

  tolerance = 1e-5;
  s = adapt_wavelet ({h}, &tolerance, 10, list = {h,u});
  fprintf (ferr, "refined: %d coarsened: %d\n", s.nf, s.nc);

  foreach()
    div[] = u.x[] - u.x[1,0] + u.y[] - u.y[0,1];
  sdiv = statsf(div);
  fprintf (ferr, "div after: %g\n", sdiv.max);

#if DEBUG
  fp = fopen ("cells1", "w");
  output_cells (fp);
  fclose (fp);
#endif

  if (error())
    return 1;

  free_grid();
}
