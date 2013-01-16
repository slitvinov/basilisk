#include "grid.h"
#include "utils.c"
#include <time.h>
#include "init.h"
#include <assert.h>

void restriction (Data * m, int n, var a)
{
  foreach_fine_to_coarse (m, n) {
    coarse(a,0,0) = (fine(a,0,0) + fine(a,1,0) + fine(a,0,1) + fine(a,1,1))/4.;
  } end_foreach_fine_to_coarse();
}

void error (Data * m, int n, var a, var e)
{
  foreach_fine_to_coarse (m, n) {
    double eb = (coarse(a,0,0) + coarse(a,0,-1))/2. - (fine(a,0,0) + fine(a,0,-1) + 
						       fine(a,1,0) + fine(a,1,-1))/4.;
    double et = (coarse(a,0,0) + coarse(a,0, 1))/2. - (fine(a,0,1) + fine(a,1,1) + 
						       fine(a,0,2) + fine(a,1,2))/4.;
    double el = (coarse(a,0,0) + coarse(a,-1,0))/2. - (fine(a,0,0) + fine(a,-1,0) + 
						       fine(a,0,1) + fine(a,-1,1))/4.;
    double er = (coarse(a,0,0) + coarse(a,1, 0))/2. - (fine(a,1,0) + fine(a,1,1) + 
						       fine(a,2,0) + fine(a,2,1))/4.;
    double avg = (eb + et + er + el)/4.;
    fine(e,0,0) = fine(e,0,1) = fine(e,1,0) = fine(e,1,1) = avg/4.;
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
