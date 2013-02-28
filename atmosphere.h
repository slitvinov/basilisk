#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "utils.h"
#include "events.h"

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
void initial_conditions (void);
void boundary_b         (void);
void boundary_h         (scalar h);
void boundary_u         (scalar u, scalar v);

void advection_centered (scalar f, scalar u, scalar v, scalar df)
{
  foreach()
    df[] = ((f[] + f[-1,0)]*u[] - 
	    (f[] + f[1,0)]*u[1,0] +
	    (f[] + f[0,-1)]*v[] - 
	    (f[] + f[0,1)]*v[0,1)]/(2.*L0*delta);
}

void advection_upwind (scalar f, scalar u, scalar v, scalar df)
{
  foreach()
    df[] = ((u[] < 0. ? f[] : f[-1,0)]*u[] - 
	    (u[1,0] > 0. ? f[] : f[1,0)]*u[1,0] +
	    (v[] < 0. ? f[] : f[0,-1)]*v[] - 
	    (v[0,1] > 0. ? f[] : f[0,1)]*v[0,1)]/(L0*delta);
}

double timestep (void)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach() {
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

void momentum (scalar u, scalar v, scalar h, scalar du, scalar dv)
{
  foreach() {
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

void ke_psi (scalar u, scalar v)
{
  foreach() {
#if 1
    ke[] = (sq(u[] + u[1,0]) + sq(v[] + v[0,1]))/8.;
#else
    double uc = u[]*u[] + u[1,0]*u[1,0];
    double vc = v[]*v[] + v[0,1]*v[0,1];
    ke[] = (uc + vc)/4.;
#endif
    psi[] = (v[] - v[-1,0] + u[0,-1] - u[])/DX;
  }
  foreach_boundary (top)
    psi[0,1] = (v[0,1] - v[-1,1] + u[] - u[0,1])/DX;
  foreach_boundary (right)
    psi[1,0] = (v[1,0] - v[] + u[1,-1] - u[1,0])/DX;
  foreach_boundary (left)
    ke[-1,0] = (sq(u[-1,0] + u[]) + sq(v[-1,0] + v[-1,1]))/8.;
  foreach_boundary (bottom)
    ke[0,-1] = (sq(u[0,-1] + u[1,-1]) + sq(v[0,-1] + v[]))/8.;
}

void advance (double t, scalar * f, scalar * df)
{
  scalar u = f[0], v = f[1], h = f[2];
  scalar du = df[0], dv = df[1], dh = df[2];

  advection_centered (h, u, v, dh);
  momentum (u, v, h, du, dv);
}

void update (double t, scalar * f)
{
  scalar u = f[0], v = f[1], h = f[2];
  boundary_h (h);
  boundary_u (u, v);
  ke_psi (u, v);
}

void run (void)
{
  parameters ();
  init_grid (N);
  events_init ();
  initial_conditions ();
  boundary_b ();
  boundary_h (h);
  boundary_u (u, v);
  ke_psi (u, v);

  clock_t start = clock ();
  struct timeval tvstart;
  gettimeofday (&tvstart, NULL);
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep ());
#if 1
    advection_centered (h, u, v, hn);
    foreach() { h[] += hn[]*dt; }
    boundary_h (h);
    momentum (u, v, h, un, vn);
    foreach() {
      u[] += un[]*dt;
      v[] += vn[]*dt;
    }
    boundary_u (u, v);
    ke_psi (u, v);
#else /* unstable! */
    scalar f[3] = { u, v, h };
    scalar df[2][3] = {{ un,  vn,  hn },
		       { un1, vn1, hn1 }};
    runge_kutta (2, t, dt, 3, f, df, advance, update);
#endif
    i++; t = tnext;
  }
  clock_t end = clock ();
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  double real = (tvend.tv_sec - tvstart.tv_sec) + (tvend.tv_usec - tvstart.tv_usec)/1e6;
  fprintf (stderr, "# " GRIDNAME ", %d steps, %g CPU, %.4g real, %.3g points.steps/s\n",
	   i, cpu, real, (N*N*(double)i/real));

  free_grid ();
}
