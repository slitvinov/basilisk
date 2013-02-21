// #include "grid/cartesian.h"
#include "advection.h"

int refine_right (Point point, void * data)
{
  VARIABLES;
  return x < -0.1;
}

int refine_circle (Point point, void * data)
{
  VARIABLES;
  return x*x + y*y < 0.25*0.25;
}

void boundary_f (void * grid, var f)
{
#if 0
  restriction (grid, f, f);
  update_halo (grid, -1, f, f);
#endif
}

void boundary_u_v (void * grid, var u, var v)
{
#if 0
  restriction_u_v (grid, u, v);
  update_halo_u_v (grid, -1, u, v);
#endif
}

void boundary_gradient (void * grid, var fx, var fy)
{
#if 0
  restriction (grid, fx, fy);
  update_halo (grid, -1, fx, fy);
#endif
}

void parameters ()
{
  // number of grid points
  N = 128;
  // end time
  TMAX = 5;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

void initial_conditions (void * grid)
{
  //  refine_function (grid, 0, -1, refine_circle, NULL);
  //  flag_halo_cells (grid);

  foreach (grid)
    f[] = exp(-100.*(sq(x + 0.2) + sq(y + .236338)));
}

int events (void * grid, int i, double t, double dt)
{
  if (i % 20 == 0)
    output_matrix (grid, f, N, stdout);

  double sum = 0., min = 1e100, max = -1e100;
  foreach (grid) {
    sum += f[]*delta*delta;
    if (f[] > max) max = f[];
    if (f[] < min) min = f[];
  }
  fprintf (stderr, "%f %f %.12f %g %g\n", t, dt, sum, min, max);

  foreach (grid) {
    u[] =   1.5*sin(2.*pi*t/5.)*sin((xu + 0.5)*pi)*cos((yu + 0.5)*pi);
    v[] = - 1.5*sin(2.*pi*t/5.)*cos((xv + 0.5)*pi)*sin((yv + 0.5)*pi);
  }
  boundary_u_v (grid, u, v);

  return 0; /* continue */
}

int main() { run (); }
