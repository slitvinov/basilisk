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
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
}

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

void initial_conditions (void * grid)
{
  //  refine_function (grid, 0, -1, refine_circle, NULL);
  //  flag_halo_cells (grid);

  foreach (grid)
    f[] = bump(x,y);
}

event (t += 0.1; t <= 5.) output_matrix (grid, f, N, stdout);

event (i++) {
  double sum = 0., min = 1e100, max = -1e100;
  foreach (grid) {
    sum += f[]*delta*delta;
    if (f[] > max) max = f[];
    if (f[] < min) min = f[];
  }
  fprintf (stderr, "%f %.12f %g %g\n", t, sum, min, max);

  foreach (grid) {
    double u1 = 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*cos((y + 0.5)*pi);
    double u0 = 1.5*sin(2.*pi*t/5.)*sin(((x - delta) + 0.5)*pi)*cos((y + 0.5)*pi);
    u[] = (u0 + u1)/2.; 
    double v1 = - 1.5*sin(2.*pi*t/5.)*cos((x + 0.5)*pi)*sin((y + 0.5)*pi);
    double v0 = - 1.5*sin(2.*pi*t/5.)*cos((x + 0.5)*pi)*sin(((y - delta) + 0.5)*pi);
    v[] = (v0 + v1)/2.;
  }
  boundary_u_v (grid, u, v);
}

event (t = 5) {
  double max = 0., norm1 = 0., norm2 = 0., area = 0.;
  var e = new var;
  foreach (grid) {
    e[] = f[] - bump(x,y);
    if (fabs(e[]) > max) max = fabs(e[]);
    double a = sq(delta);
    norm1 += a*fabs(e[]);
    norm2 += a*e[]*e[];
    area  += a;
  }
  fprintf (stderr, "error: %g %g %g\n", norm1/area, sqrt(norm2/area), max);

  FILE * fp = fopen ("error.m", "w");
  output_matrix (grid, e, N, fp);
  fclose (fp);
}

int main() { run (); }
