#include <time.h>
#include "utils.h"
#include "events.h"

scalar f = new scalar, u = new scalar, v = new scalar;

// user-provided functions
void parameters         (void);
void initial_conditions (void * grid);
void boundary_f         (void * grid, scalar f);
void boundary_u_v       (void * grid, scalar u, scalar v);
void boundary_gradient  (void * grid, scalar fx, scalar fy);

double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod   (double r) { return generic_limiter (r, 1.);  }
double superbee (double r) { return generic_limiter (r, 2.);  }
double sweby    (double r) { return generic_limiter (r, 1.5); }

void gradient (void * grid, const scalar f, vector g)
{
  foreach (grid)
    foreach_dimension ()
      g.x[] = (f[1,0] - f[-1,0])/(2.*delta);
  //      g.x[] = minmod ((f[1,0] - f[])/(f[] - f[-1,0]))*(f[] - f[-1,0])/delta;
}

void fluxes_upwind_bcg (void * grid,
			const scalar f, const vector g,
			const vector u, 
			vector flux,
			double dt)
{
  foreach (grid)
    foreach_dimension () {
      double un = dt*u.x[]/delta, s = sign(un);
      int i = -(s + 1.)/2.;
      double f2 = f[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*delta/2.;
      double vn = u.y[i,0] + u.y[i,1];
      double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
      f2 -= dt*vn*fyy/(4.*delta);
      flux.x[] = f2*u.x[];
    }
}

double timestep (void * grid, const scalar u, const scalar v)
{
  double dtmax = DT/CFL;
  foreach (grid) {
    double dx = L0*delta;
    double dt = dx/fabs(u[]);
    if (dt < dtmax) dtmax = dt;
    dt = dx/fabs(v[]);
    if (dt < dtmax) dtmax = dt;
  }
  return dtmax*CFL;
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
    double dt = dtnext (t, timestep (grid, u, v));
    vector flux = new vector, g = new vector, uv = {u,v};
    gradient (grid, f, g);
    boundary_gradient (grid, g.x, g.y);
    fluxes_upwind_bcg (grid, f, g, uv, flux, dt);
    boundary_u_v (grid, flux.x, flux.y);
    foreach (grid)
      f[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/delta;
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
