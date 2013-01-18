#include <time.h>
#include <assert.h>

struct _Data {
  double h, u, v, e, w;
};
#define h(k,l) data(k,l).h
#define e(k,l) data(k,l).e

#include "utils.h"
#include "quadtree.c"
#include "utils.c"

#include "init.h"

void restriction (Data * m, int n, var a)
{
  foreach_fine_to_coarse (m, n) {
    stencil(a,0,0) = (fine(a,0,0) + fine(a,1,0) + fine(a,0,1) + fine(a,1,1))/4.;
  } end_foreach_fine_to_coarse();
}

void wavelet (Data * m, int n, var a, var w)
{
  foreach_fine_to_coarse (m, n) {
    /* difference between fine value and bilinearly-interpolated coarse value */
    fine(w,0,0) = fine(a,0,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,-1)) + coarse(a,-1,-1))/16.;
    fine(w,0,1) = fine(a,0,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,1)) + coarse(a,-1,1))/16.;
    fine(w,1,0) = fine(a,1,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,0,-1) + coarse(a,1,0)) + coarse(a,1,-1))/16.;
    fine(w,1,1) = fine(a,1,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,1,0) + coarse(a,0,1)) + coarse(a,1,1))/16.;
  } end_foreach_fine_to_coarse ();
}

void coarsen_wavelet (Data * m, int n, var w, double max)
{
  foreach_fine_to_coarse (m, n) {
    double error = 0.;
    foreach_fine(k,l) {
      double e = fabs(fine(w,k,l));
      if (e > error)
	error = e;
    }
    if (error < max) {
      /* coarsen */
      cell.flags |= leaf;
      foreach_fine(k,l) {
	finecell(k,l).flags &= !leaf;
	finecell(k,l).flags |= inactive;
      }
    }
    /* propagate the error to coarser levels */
    coarse(w,0,0) = fabs(coarse(w,0,0)) + error;
  } end_foreach_fine_to_coarse ();
}

void coarsen (Data * m, int n, var cost, double max)
{
  foreach_cell (m, n)
    if (cell.flags & leaf)
      continue;
    else if (val(cost) < max) {
      /* coarsen */
      cell.flags |= leaf;
      continue;
    }
  end_foreach_cell();
}

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = atoi(argv[1]);
  void * m = init_grid (n);

  initial_conditions (m, n);
  symmetry (m, n, var(h));
  uv_symmetry (m, n, var(u), var(v));
  
  /* initial coarsening */
  restriction (m, n, var(h));
  wavelet (m, n, var(h), var(w));
  coarsen_wavelet (m, n, var(w), 1e-4);
  //  coarsen (m, n, var(e), 1e-4);

  foreach_cell (m , n)
    if (cell.flags & inactive)
      printf ("%g %g %d inactive\n", x, y, level);
    else if (cell.flags & leaf)
      printf ("%g %g %d leaf\n", x, y, level);
  end_foreach_cell();

  clock_t start, end;
  start = clock ();
  do {
    //    double dt = timestep (m, n);
    double dt = CFL/n;
    //    #include "output.h"
    //    tracer_advection_upwind (m, n, dt);
    //    boundary_h (m, n);
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRID ", %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));

  free_grid (m);
}
