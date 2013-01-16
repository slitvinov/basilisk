#include "grid.h"
#include "utils.c"
#include <time.h>
#include "init.h"
#include <assert.h>

void restriction (Data * m, int n)
{
  foreach_fine_to_coarse (m, n) {
    h(0,0) = (fine(0,0).h + fine(1,0).h + fine(0,1).h + fine(1,1).h)/4.;
  } end_foreach_fine_to_coarse();
}

void error (Data * m, int n)
{
  foreach_fine_to_coarse (m, n) {
    double eb = (h(0,0) + h(0,-1))/2. - (fine(0,0).h + fine(0,-1).h + 
					 fine(1,0).h + fine(1,-1).h)/4.;
    double et = (h(0,0) + h(0, 1))/2. - (fine(0,1).h + fine(1,1).h + 
					 fine(0,2).h + fine(1,2).h)/4.;
    double el = (h(0,0) + h(-1,0))/2. - (fine(0,0).h + fine(-1,0).h + 
					 fine(0,1).h + fine(-1,1).h)/4.;
    double er = (h(0,0) + h(1, 0))/2. - (fine(1,0).h + fine(1,1).h + 
					 fine(2,0).h + fine(2,1).h)/4.;
#if 0
    if (l == 7)
      printf ("%g %g %g %g %g %g %g\n", XC(i), YC(j), eb, et, er, el,
	      (eb + et + er + el)/4.);
#endif
    fine(0,0).b = fine(0,1).b = fine(1,0).b = fine(1,1).b = (eb + et + er + el)/4.;
  } end_foreach_fine_to_coarse();
}

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = N;

  Data * m = init_grid (n);
  initial_conditions (m, n);
  boundary_h (m, n);
  boundary_u (m, n);

  clock_t start, end;
  start = clock ();
  do {
    double dt = timestep (m, n);
    //    double dt = CFL/n;
    #include "output.h"
    tracer_advection_upwind (m, n, dt);
    boundary_h (m, n);
#if 0
    restriction (m, n);
    error (m, n);
#endif
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
