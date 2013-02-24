#include <time.h>
#include "utils.h"
#include "events.h"

var f = new var, u = new var, v = new var;

// user-provided functions
void parameters         (void);
void initial_conditions (void * grid);
void boundary_f         (void * grid, var f);
void boundary_u_v       (void * grid, var u, var v);
void boundary_gradient  (void * grid, var fx, var fy);

void fluxes_upwind (void * grid, const var f, const var u, const var v, var fu, var fv)
{
  foreach (grid) {
    fu[] = (u[] < 0. ? f[] : f[-1,0])*u[];
    fv[] = (v[] < 0. ? f[] : f[0,-1])*v[];
  }
}

void gradient (void * grid, const var f, var fx, var fy)
{
  foreach (grid) {
    delta *= 2.;
    fx[] = (f[1,0] - f[-1,0])/delta;
    fy[] = (f[0,1] - f[0,-1])/delta;
  }
}

double limiter (double r)
{
  return 1.; // no limiter
  return r*(r + 1.)/(r*r + 1.); // van Albada
  
  //  double beta = 1.;
  //  return max(0, max(min(beta*r, 1.), min(r, beta))); // Sweby
  return max(0., min(1., r)); // minmod
}

void fluxes_upwind_bcg (void * grid,
			const var f, const var fx, const var fy,
			const var u, const var v,
			var fu, var fv,
			double dt)
{
  foreach (grid) {
    double r, f1, f2;
    if (u[] < 0.) {
      f1 = f[];
      double un = dt*(u[] + u[1,0])/(2.*delta);
      f2 = f1 + max(-1., -1. - un)*fx[]*delta/2.;
      double vn = v[] + v[0,1];
      double fyy = vn > 0. ? f[0,1] - f[] : f[] - f[0,-1];
      f2 -= dt*vn*fyy/(4.*delta);
    }
    else {
      f1 = f[-1,0];
      double un = dt*(u[-1,0] + u[])/(2.*delta);
      f2 = f1 + min(1., 1. - un)*fx[-1,0]*delta/2.;
      double vn = v[-1,0] + v[-1,1];
      double fyy = vn > 0. ? f[-1,1] - f[-1,0] : f[-1,0] - f[-1,-1];
      f2 -= dt*vn*fyy/(4.*delta);
    }
    r = (f[] - f[-1,0])/(f[1,0] - f[]);
    fu[] = (f1 - limiter(r)*(f1 - f2))*u[];
    
    if (v[] < 0.) {
      f1 = f[];
      double un = dt*(v[] + v[0,1])/(2.*delta);
      f2 = f1 + max(-1., -1. - un)*fy[]*delta/2.;
      double vn = u[] + u[1,0];
      double fxx = vn > 0. ? f[1,0] - f[] : f[] - f[-1,0];
      f2 -= dt*vn*fxx/(4.*delta);
    }
    else {
      f1 = f[0,-1];
      double un = dt*(v[0,-1] + v[])/(2.*delta);
      f2 = f1 + min(1., 1. - un)*fy[0,-1]*delta/2.;
      double vn = u[0,-1] + u[1,-1];
      double fxx = vn > 0. ? f[1,-1] - f[0,-1] : f[0,-1] - f[-1,-1];
      f2 -= dt*vn*fxx/(4.*delta);
    }
    r = (f[] - f[0,-1])/(f[0,1] - f[]);
    fv[] = (f1 - limiter(r)*(f1 - f2))*v[];    
  }
}

double timestep (void * grid, double t, const var u, const var v)
{
  double dtmax = DT/CFL;
  foreach (grid) {
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
  dtmax *= CFL;
  int n = (tnext - t)/dtmax;
  if (n == 0)
    dtmax = tnext - t;
  else {
    dtmax = (tnext - t)/(n + 1);
    tnext = t + dtmax;
  }
  return dtmax;
}

void run (void)
{
  parameters ();
  void * grid = init_grid (N);

  events_init (grid);
  initial_conditions (grid);
  boundary_f (grid, f);

  clock_t cstart = clock ();
  double t = 0.;
  int i = 0;
  while (events (grid, i, t)) {
    double dt = timestep (grid, t, u, v);
    var fu = new var, fv = new var;
    var fx = new var, fy = new var;
    gradient (grid, f, fx, fy);
    boundary_gradient (grid, fx, fy);
    fluxes_upwind_bcg (grid, f, fx, fy, u, v, fu, fv, dt);
    boundary_u_v (grid, fu, fv);
    foreach (grid)
      f[] += dt*(fu[] - fu[1,0] + fv[] - fv[0,1])/delta;
    boundary_f (grid, f);
    i++; t = tnext;
  }
  clock_t cend = clock ();
  double cpu = ((double) (cend - cstart))/CLOCKS_PER_SEC;
  int n = 0;
  foreach (grid) n++;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (n*(double)i/cpu));

  free_grid (grid);
}
