#include <time.h>
#include "utils.h"

new var u, v, h, b, ke, psi, un, vn, hn;

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;
// Viscosity
double NU = 0.;

void advection_centered (void * grid, var f, var u, var v, var df)
{
  foreach (grid)
    df(0,0) = ((f(0,0) + f(-1,0))*u(0,0) - 
	       (f(0,0) + f(1,0))*u(1,0) +
	       (f(0,0) + f(0,-1))*v(0,0) - 
	       (f(0,0) + f(0,1))*v(0,1))/(2.*L0*delta);
}

void advection_upwind (void * grid, var f, var u, var v, var df)
{
  foreach (grid)
    df(0,0) = ((u(0,0) < 0. ? f(0,0) : f(-1,0))*u(0,0) - 
	       (u(1,0) > 0. ? f(0,0) : f(1,0))*u(1,0) +
	       (v(0,0) < 0. ? f(0,0) : f(0,-1))*v(0,0) - 
	       (v(0,1) > 0. ? f(0,0) : f(0,1))*v(0,1))/(L0*delta);
}

double timestep (void * grid)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach (grid) {
    double dx = L0*delta;
    dx *= dx;
    if (h(0,0) > 0.) {
      double dt = dx/(G*h(0,0));
      if (dt < dtmax) dtmax = dt;
    }
    if (u(0,0) != 0.) {
      double dt = dx/(u(0,0)*u(0,0));
      if (dt < dtmax) dtmax = dt;
    }
    if (v(0,0) != 0.) {
      double dt = dx/(v(0,0)*v(0,0));
      if (dt < dtmax) dtmax = dt;
    }
  }
  return sqrt (dtmax)*CFL;
}

void momentum (void * grid, var u, var v, var h, var du, var dv)
{
  foreach (grid) {
    double g = G*(h(0,0) + b(0,0)) + ke(0,0);
    double psiu = (psi(0,0) + psi(0,1))/2.;
    double dx = L0*delta;
    du(0,0) = 
      - (g - G*(h(-1,0) + b(-1,0)) - ke(-1,0))/dx
      + (psiu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1))/4.
      + NU*(u(1,0) + u(0,1) + u(-1,0) + u(0,-1) - 4.*u(0,0))/(dx*dx);
    double psiv = (psi(0,0) + psi(1,0))/2.;
    dv(0,0) = 
      - (g - G*(h(0,-1) + b(0,-1)) - ke(0,-1))/dx
      - (psiv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1))/4.
      + NU*(v(1,0) + v(0,1) + v(-1,0) + v(0,-1) - 4.*v(0,0))/(dx*dx);
  }
}

void ke_psi (void * grid, var u, var v)
{
  foreach (grid) {
#if 1
    ke(0,0) = (sq(u(0,0) + u(1,0)) + sq(v(0,0) + v(0,1)))/8.;
#else
    double uc = u(0,0)*u(0,0) + u(1,0)*u(1,0);
    double vc = v(0,0)*v(0,0) + v(0,1)*v(0,1);
    ke(0,0) = (uc + vc)/4.;
#endif
    psi(0,0) = (v(0,0) - v(-1,0) + u(0,-1) - u(0,0))/DX;
  }
  foreach_boundary (grid, top)
    psi(0,1) = (v(0,1) - v(-1,1) + u(0,0) - u(0,1))/DX;
  foreach_boundary (grid, right)
    psi(1,0) = (v(1,0) - v(0,0) + u(1,-1) - u(1,0))/DX;
  foreach_boundary (grid, left)
    ke(-1,0) = (sq(u(-1,0) + u(0,0)) + sq(v(-1,0) + v(-1,1)))/8.;
  foreach_boundary (grid, bottom)
    ke(0,-1) = (sq(u(0,-1) + u(1,-1)) + sq(v(0,-1) + v(0,0)))/8.;
}

void advance (void * grid, double t, var * f, var * df)
{
  var u = f[0], v = f[1], h = f[2];
  var du = df[0], dv = df[1], dh = df[2];

  advection_centered (grid, h, u, v, dh);
  momentum (grid, u, v, h, du, dv);
}

#include "parameters.h"

void update (void * grid, double t, var * f)
{
  var u = f[0], v = f[1], h = f[2];
  boundary_h (grid, h);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);
}

int main (int argc, char ** argv)
{
  parameters();
  double t = 0;
  int i = 0;

  void * grid = init_grid (N);
  initial_conditions (grid);
  boundary_b (grid);
  boundary_h (grid, h);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);

  clock_t start, end;
  start = clock ();
  do {
    double dt = timestep (grid);
    events(grid, i, t, dt);
#if 1
    advection_centered (grid, h, u, v, hn);
    foreach (grid) { h(0,0) += hn(0,0)*dt; }
    boundary_h (grid, h);
    momentum (grid, u, v, h, un, vn);
    foreach (grid) {
      u(0,0) += un(0,0)*dt;
      v(0,0) += vn(0,0)*dt;
    }
    boundary_u (grid, u, v);
    ke_psi (grid, u, v);
#else /* unstable! */
    var f[3] =  { u, v, h };
    var df[2][3] = {{ un,  vn,  hn },
		    { un1, vn1, hn1 }};
    runge_kutta (2, grid, t, dt, 3, f, df, advance, update);
#endif
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (N*N*(double)i/cpu));

  free_grid (grid);
}
