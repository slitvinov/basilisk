/* tangential interpolation on staggered vector fields  */

#include "utils.h"

scalar h[];
staggered vector u[];

static double solution (double x, double y)
{
#if 0
  double R0 = 0.1;
  return exp(-(x*x+y*y)/sq(R0));
#else
  return x + y;
#endif
}

u.x[right]  = solution(x,y);
u.x[left]   = solution(x,y);
u.y[top]    = solution(x,y);
u.y[bottom] = solution(x,y);

u.x[top]    = dirichlet(solution(x,y));
u.x[bottom] = dirichlet(solution(x,y));
u.y[right]  = dirichlet(solution(x,y));
u.y[left]   = dirichlet(solution(x,y));

double tolerance = 1e-4;

static int error()
{
  double max = 0, maxv = 0.;

  scalar eu[];
  foreach()
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

  fprintf (stderr, "maximum error: %g %g\n", max, maxv);
  stats s = statsf (eu);
  fprintf (stderr, "eu: avg: %g stddev: %g max: %g\n", 
	   s.sum/s.area, s.stddev, s.max);

  return (max != maxv);
}

int main (int argc, char ** argv)
{
  int n = 1024;
  init_grid (n);

  X0 = Y0 = -0.5;
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/sq(R0));
  boundary ({h});

  trash ({u});
  foreach_face()
    u.x[] = solution(x,y);
  boundary ((scalar *){u});

  astats s = adapt_wavelet ({h}, &tolerance, 10);
  fprintf (stderr, "refined: %d coarsened: %d\n", s.nf, s.nc);

#if DEBUG
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);
#endif

  if (error())
    return 1;

  tolerance = 1e-5;
  s = adapt_wavelet ({h}, &tolerance, 10);
  fprintf (stderr, "refined: %d coarsened: %d\n", s.nf, s.nc);

#if DEBUG
  fp = fopen ("cells1", "w");
  output_cells (fp);
  fclose (fp);
#endif

  if (error())
    return 1;

  free_grid();
}