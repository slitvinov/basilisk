#include <time.h>
#include "utils.h"
#include "mg.h"

scalar u = new scalar, v = new scalar, p = new scalar;
scalar Sxx = new scalar, Syy = new scalar, Sxy = new scalar;

// Default parameters, do not change them!! edit parameters.h instead
// Viscosity
double NU = 0.;
// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on maximum divergence
double TOLERANCE = 1e-3;
// user-provided functions
void parameters         (void);
void initial_conditions (void);
void boundary_p         (scalar p, int l);
void boundary_u         (scalar u, scalar v);
int  events             (int i, double t, double dt);
void end                (void);

double timestep ()
{
  double dtmax = DT/CFL;
  foreach() {
    double dx = L0*delta;
    if (u[] != 0.) {
      double dt = dx/fabs(u[]);
      if (dt < dtmax) dtmax = dt;
    }
    if (v[] != 0.) {
      double dt = dx/fabs(v[]);
      if (dt < dtmax) dtmax = dt;
    }
  }
  return dtmax*CFL;
}

#define sq(x) ((x)*(x))

void stresses ()
{
  foreach() {
    delta *= L0;
    Sxx[] = 
      - sq(u[] + u[1,0])/4. 
      + 2.*NU*(u[1,0] - u[])/delta;
    Syy[] = 
      - sq(v[] + v[0,1])/4. 
      + 2.*NU*(v[0,1] - v[])/delta;
    Sxy[] = 
      - (u[] + u[0,-1])*(v[] + v[-1,0])/4. +
      NU*(u[] - u[0,-1] + v[] - v[-1,0])/delta;
  }
  foreach_boundary (left) {
    delta *= L0;
    Sxx[-1,0] = 
      - sq(u[-1,0] + u[])/4. 
      + 2.*NU*(u[] - u[-1,0])/delta;
  }
  foreach_boundary (top) {
    delta *= L0;
    Sxy[0,1] = 
      - (u[0,1] + u[])*(v[0,1] + v[-1,1])/4. +
      NU*(u[0,1] - u[] + v[0,1] - v[-1,1])/delta;
  }
  foreach_boundary (right) {
    delta *= L0;
    Sxy[1,0] = 
      - (u[1,0] + u[1,-1])*(v[1,0] + v[])/4. +
      NU*(u[1,0] - u[1,-1] + v[1,0] - v[])/delta;
  }
  foreach_boundary (bottom) {
    delta *= L0;
    Syy[0,-1] = 
      - sq(v[0,-1] + v[])/4. 
      + 2.*NU*(v[] - v[0,-1])/delta;
  }
}

void advance (double dt)
{
  foreach() {
    delta *= L0;
    u[] += dt*(Sxx[] - Sxx[-1,0] + Sxy[0,1] - Sxy[])/delta;
    v[] += dt*(Syy[] - Syy[0,-1] + Sxy[1,0] - Sxy[])/delta;
  }
}

void relax (scalar a, scalar b, int l)
{
  foreach_level (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - L0*L0*delta*delta*b[])/4.;
}

double residual (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
  foreach() {
    res[] = b[] + (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/(L0*L0*delta*delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

void projection (scalar u, scalar v, scalar p, 
		 scalar div, scalar res, scalar dp)
{
  double sum = 0.;
  foreach() {
    div[] = (u[1,0] - u[] + v[0,1] - v[])/(L0*delta);
    sum += div[];
  }
  double maxres = residual (p, div, res);
  int i;
  for (i = 0; i < NITERMAX && (i < 1 || maxres > TOLERANCE); i++) {
    mg_cycle (p, res, dp,
	      relax, boundary_p,
	      4, 0);
    boundary_p (p, depth());
    maxres = residual (p, div, res);
  }
  if (i == NITERMAX)
    fprintf (stderr, "WARNING: convergence not reached after %d iterations\n  sum: %g\n", 
	     NITERMAX, sum);
  foreach() {
    delta *= L0;
    u[] -= (p[] - p[-1,0])/delta;
    v[] -= (p[] - p[0,-1])/delta;
  }
}

void run (void)
{
  parameters ();

  double t = 0;
  int i = 0;

  init_grid (N);
  initial_conditions ();
  boundary_u (u, v);
  projection (u, v, p, Sxx, Syy, Sxy);
  boundary_u (u, v);

  clock_t cstart = clock ();
  do {
    double dt = timestep ();
    if (events (i, t, dt))
      break;
    stresses ();
    advance (dt);
    boundary_u (u, v);
    projection (u, v, p, Sxx, Syy, Sxy);
    boundary_u (u, v);
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end ();
  clock_t cend = clock ();
  double cpu = ((double) (cend - cstart))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (N*N*(double)i/cpu));

  free_grid ();
}
