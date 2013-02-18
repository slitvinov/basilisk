#include <time.h>
#include "utils.h"
#include "mg.h"

var u = new var, v = new var, p = new var;
var Sxx = new var, Syy = new var, Sxy = new var;

// Default parameters, do not change them!! edit parameters.h instead
// Viscosity
double NU = 0.;
// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on maximum divergence
double TOLERANCE = 1e-3;
// user-provided functions
void parameters         (void);
void initial_conditions (void * grid);
void boundary_p         (void * grid, var p, int l);
void boundary_u         (void * grid, var u, var v);
int  events             (void * grid, int i, double t, double dt);
void end                (void * grid);

double timestep (void * grid)
{
  double dtmax = DT/CFL;
  foreach (grid) {
    double dx = L0*delta;
    if (u(0,0) != 0.) {
      double dt = dx/fabs(u(0,0));
      if (dt < dtmax) dtmax = dt;
    }
    if (v(0,0) != 0.) {
      double dt = dx/fabs(v(0,0));
      if (dt < dtmax) dtmax = dt;
    }
  }
  return dtmax*CFL;
}

#define sq(x) ((x)*(x))

void stresses (void * grid)
{
  foreach (grid) {
    delta *= L0;
    Sxx(0,0) = 
      - sq(u(0,0) + u(1,0))/4. 
      + 2.*NU*(u(1,0) - u(0,0))/delta;
    Syy(0,0) = 
      - sq(v(0,0) + v(0,1))/4. 
      + 2.*NU*(v(0,1) - v(0,0))/delta;
    Sxy(0,0) = 
      - (u(0,0) + u(0,-1))*(v(0,0) + v(-1,0))/4. +
      NU*(u(0,0) - u(0,-1) + v(0,0) - v(-1,0))/delta;
  }
  foreach_boundary (grid, left) {
    delta *= L0;
    Sxx(-1,0) = 
      - sq(u(-1,0) + u(0,0))/4. 
      + 2.*NU*(u(0,0) - u(-1,0))/delta;
  }
  foreach_boundary (grid, top) {
    delta *= L0;
    Sxy(0,1) = 
      - (u(0,1) + u(0,0))*(v(0,1) + v(-1,1))/4. +
      NU*(u(0,1) - u(0,0) + v(0,1) - v(-1,1))/delta;
  }
  foreach_boundary (grid, right) {
    delta *= L0;
    Sxy(1,0) = 
      - (u(1,0) + u(1,-1))*(v(1,0) + v(0,0))/4. +
      NU*(u(1,0) - u(1,-1) + v(1,0) - v(0,0))/delta;
  }
  foreach_boundary (grid, bottom) {
    delta *= L0;
    Syy(0,-1) = 
      - sq(v(0,-1) + v(0,0))/4. 
      + 2.*NU*(v(0,0) - v(0,-1))/delta;
  }
}

void advance (void * grid, double dt)
{
  foreach (grid) {
    delta *= L0;
    u(0,0) += dt*(Sxx(0,0) - Sxx(-1,0) + Sxy(0,1) - Sxy(0,0))/delta;
    v(0,0) += dt*(Syy(0,0) - Syy(0,-1) + Sxy(1,0) - Sxy(0,0))/delta;
  }
}

void relax (void * grid, var a, var b, int l)
{
  foreach_level (grid, l)
    a(0,0) = (a(1,0) + a(-1,0) +
		  a(0,1) + a(0,-1) 
		  - L0*L0*delta*delta*b(0,0))/4.;
}

double residual (void * grid, var a, var b, var res)
{
  double maxres = 0.;
  foreach (grid) {
    res(0,0) = b(0,0) + 
    (4.*a(0,0) - a(1,0) - a(-1,0) - a(0,1) - a(0,-1))/(L0*L0*delta*delta);
    if (fabs (res(0,0)) > maxres)
      maxres = fabs (res(0,0));
  }
  return maxres;
}

void projection (void * grid, var u, var v, var p, 
		 var div, var res, var dp)
{
  double sum = 0.;
  foreach(grid) {
    div(0,0) = (u(1,0) - u(0,0) + v(0,1) - v(0,0))/(L0*delta);
    sum += div(0,0);
  }
  double maxres = residual (grid, p, div, res);
  int i;
  for (i = 0; i < NITERMAX && (i < 1 || maxres > TOLERANCE); i++) {
    mg_cycle (grid, depth(grid), p, res, dp,
	      relax, boundary_p,
	      4, 0);
    boundary_p (grid, p, depth(grid));
    maxres = residual (grid, p, div, res);
  }
  if (i == NITERMAX)
    fprintf (stderr, "WARNING: convergence not reached after %d iterations\n  sum: %g\n", 
	     NITERMAX, sum);
  foreach (grid) {
    delta *= L0;
    u(0,0) -= (p(0,0) - p(-1,0))/delta;
    v(0,0) -= (p(0,0) - p(0,-1))/delta;
  }
}

void run (void)
{
  parameters ();

  double t = 0;
  int i = 0;

  void * grid = init_grid (N);
  initial_conditions (grid);
  boundary_u (grid, u, v);
  projection (grid, u, v, p, Sxx, Syy, Sxy);
  boundary_u (grid, u, v);

  clock_t cstart, cend;
  cstart = clock ();
  do {
    double dt = timestep (grid);
    if (events (grid, i, t, dt))
      break;
    stresses (grid);
    advance (grid, dt);
    boundary_u (grid, u, v);
    projection (grid, u, v, p, Sxx, Syy, Sxy);
    boundary_u (grid, u, v);
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end (grid);
  cend = clock ();
  double cpu = ((double) (cend - cstart))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (N*N*(double)i/cpu));

  free_grid (grid);
}
