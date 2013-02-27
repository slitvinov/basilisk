#include <stdio.h>
#include <time.h>
#include "utils.h"

scalar u = new scalar, v = new scalar, h = new scalar, b = new scalar;
scalar ke = new scalar, psi = new scalar;
scalar un = new scalar, vn = new scalar, hn = new scalar;

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;
// Viscosity
double NU = 0.;
// user-provided functions
void parameters         (void);
void initial_conditions (void * grid);
void boundary_b         (void * grid);
void boundary_h         (void * grid, scalar h);
void boundary_u         (void * grid, scalar u, scalar v);
int  events             (void * grid, int i, double t, double dt);

void advection_centered (void * grid, scalar f, scalar u, scalar v, scalar df)
{
  foreach (grid)
    df[] = ((f[] + f[-1,0)]*u[] - 
	    (f[] + f[1,0)]*u[1,0] +
	    (f[] + f[0,-1)]*v[] - 
	    (f[] + f[0,1)]*v[0,1)]/(2.*L0*delta);
}

void advection_upwind (void * grid, scalar f, scalar u, scalar v, scalar df)
{
  foreach (grid)
    df[] = ((u[] < 0. ? f[] : f[-1,0)]*u[] - 
	       (u[1,0] > 0. ? f[] : f[1,0)]*u[1,0] +
	       (v[] < 0. ? f[] : f[0,-1)]*v[] - 
	       (v[0,1] > 0. ? f[] : f[0,1)]*v[0,1)]/(L0*delta);
}

double timestep (void * grid)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach (grid) {
    double dx = L0*delta;
    dx *= dx;
    if (h[] > 0.) {
      double dt = dx/(G*h[]);
      if (dt < dtmax) dtmax = dt;
    }
    if (u[] != 0.) {
      double dt = dx/(u[]*u[]);
      if (dt < dtmax) dtmax = dt;
    }
    if (v[] != 0.) {
      double dt = dx/(v[]*v[]);
      if (dt < dtmax) dtmax = dt;
    }
  }
  return sqrt (dtmax)*CFL;
}

void momentum (void * grid, scalar u, scalar v, scalar h, scalar du, scalar dv)
{
  foreach (grid) {
    double g = G*(h[] + b[]) + ke[];
    double psiu = (psi[] + psi[0,1])/2.;
    double dx = L0*delta;
    du[] = 
      - (g - G*(h[-1,0] + b[-1,0]) - ke[-1,0])/dx
      + (psiu + F0)*(v[] + v[0,1] + v[-1,0] + v[-1,1])/4.
      + NU*(u[1,0] + u[0,1] + u[-1,0] + u[0,-1] - 4.*u[])/(dx*dx);
    double psiv = (psi[] + psi[1,0])/2.;
    dv[] = 
      - (g - G*(h[0,-1] + b[0,-1]) - ke[0,-1])/dx
      - (psiv + F0)*(u[] + u[1,0] + u[0,-1] + u[1,-1])/4.
      + NU*(v[1,0] + v[0,1] + v[-1,0] + v[0,-1] - 4.*v[])/(dx*dx);
  }
}

void ke_psi (void * grid, scalar u, scalar v)
{
  foreach (grid) {
#if 1
    ke[] = (sq(u[] + u[1,0]) + sq(v[] + v[0,1]))/8.;
#else
    double uc = u[]*u[] + u[1,0]*u[1,0];
    double vc = v[]*v[] + v[0,1]*v[0,1];
    ke[] = (uc + vc)/4.;
#endif
    psi[] = (v[] - v[-1,0] + u[0,-1] - u[])/DX;
  }
  foreach_boundary (grid, top)
    psi[0,1] = (v[0,1] - v[-1,1] + u[] - u[0,1])/DX;
  foreach_boundary (grid, right)
    psi[1,0] = (v[1,0] - v[] + u[1,-1] - u[1,0])/DX;
  foreach_boundary (grid, left)
    ke[-1,0] = (sq(u[-1,0] + u[]) + sq(v[-1,0] + v[-1,1]))/8.;
  foreach_boundary (grid, bottom)
    ke[0,-1] = (sq(u[0,-1] + u[1,-1]) + sq(v[0,-1] + v[]))/8.;
}

void advance (void * grid, double t, scalar * f, scalar * df)
{
  scalar u = f[0], v = f[1], h = f[2];
  scalar du = df[0], dv = df[1], dh = df[2];

  advection_centered (grid, h, u, v, dh);
  momentum (grid, u, v, h, du, dv);
}

void update (void * grid, double t, scalar * f)
{
  scalar u = f[0], v = f[1], h = f[2];
  boundary_h (grid, h);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);
}

void run (void)
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
    foreach (grid) { h[] += hn[]*dt; }
    boundary_h (grid, h);
    momentum (grid, u, v, h, un, vn);
    foreach (grid) {
      u[] += un[]*dt;
      v[] += vn[]*dt;
    }
    boundary_u (grid, u, v);
    ke_psi (grid, u, v);
#else /* unstable! */
    scalar f[3] = { u, v, h };
    scalar df[2][3] = {{ un,  vn,  hn },
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
