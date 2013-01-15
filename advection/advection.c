#include "grid.h"
#include "utils.c"

#include <time.h>

#include "init.h"

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
    #include "output.h"
    tracer_advection_upwind (m, n, dt);
    boundary_h (m, n);
    t += dt; i++;
  } while (t < TMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
