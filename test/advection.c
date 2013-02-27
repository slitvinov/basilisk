#include "grid/cartesian.h"
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

void boundary_f (scalar f)
{
#if 0
  restriction (f, f);
  update_halo (-1, f, f);
#endif
}

void boundary_u_v (scalar u, scalar v)
{
#if 0
  restriction_u_v (u, v);
  update_halo_u_v (-1, u, v);
#endif
}

void boundary_gradient (scalar fx, scalar fy)
{
#if 0
  restriction (fx, fy);
  update_halo (-1, fx, fy);
#endif
}

void parameters ()
{
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
}

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

void initial_conditions ()
{
  //  refine_function (0, -1, refine_circle, NULL);
  //  flag_halo_cells ();

  foreach()
    f[] = bump(x,y);
}

// event (t += 0.1; t <= 5.) output_matrix (f, N, stdout);

event (i++) {
  foreach() {
    u[] = 1.5*sin(2.*pi*t/5.)*sin((xu + 0.5)*pi)*cos((yu + 0.5)*pi);
    v[] = - 1.5*sin(2.*pi*t/5.)*cos((xv + 0.5)*pi)*sin((yv + 0.5)*pi);
  }
  boundary_u_v (u, v);
}

event (t = {0,5}) {
  double sum = 0., min = 1e100, max = -1e100;
  foreach() {
    sum += f[]*delta*delta;
    if (f[] > max) max = f[];
    if (f[] < min) min = f[];
  }
  fprintf (stderr, "# %f %.12f %g %g\n", t, sum, min, max);  
}

event (t = 5) {
  double max = 0., norm1 = 0., norm2 = 0., area = 0.;
  scalar e = new scalar;
  foreach() {
    e[] = f[] - bump(x,y);
    if (fabs(e[]) > max) max = fabs(e[]);
    double a = sq(delta);
    norm1 += a*fabs(e[]);
    norm2 += a*e[]*e[];
    area  += a;
  }
  fprintf (stderr, "%d %g %g %g\n", N, norm1/area, sqrt(norm2/area), max);
  
  if (N == 256)
    output_matrix (e, N, stdout);
}

int main() {
  for (N = 64; N <= 256; N *= 2)
    run ();
}
