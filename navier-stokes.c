#include <time.h>

struct _Data {
  double u, v, p, ke, psi;
  double un, vn;
};

#include "grid/multigrid.h"
#include "utils.h"
#include "mg.h"

var u = var(u), v = var(v), p = var(p), ke = var(ke), psi = var(psi),
  un = var(un), vn = var(vn);

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 0.;
// Viscosity
double NU = 0.;
// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on maximum divergence
double TOLERANCE = 1e-3;

#include "lid.h"

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

void momentum (void * grid, var u, var v, var du, var dv)
{
  foreach (grid) {
    double g = ke(0,0);
    double psiu = (psi(0,0) + psi(0,1))/2.;
    double dx = L0*delta;
    du(0,0) = 
      - (g - ke(-1,0))/dx
      + (psiu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1))/4.
      + NU*(u(1,0) + u(0,1) + u(-1,0) + u(0,-1) - 4.*u(0,0))/(dx*dx);
    double psiv = (psi(0,0) + psi(1,0))/2.;
    dv(0,0) = 
      - (g - ke(0,-1))/dx
      - (psiv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1))/4.
      + NU*(v(1,0) + v(0,1) + v(-1,0) + v(0,-1) - 4.*v(0,0))/(dx*dx);
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

#define sq(x) ((x)*(x))

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

int main (int argc, char ** argv)
{
  parameters ();

  double t = 0;
  int i = 0;

  void * grid = init_grid (N);
  initial_conditions (grid);
  boundary_u (grid, u, v);
  projection (grid, u, v, p, un, vn, ke);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);

  clock_t cstart, cend;
  cstart = clock ();
  do {
    double dt = timestep (grid);
    if (events (grid, i, t, dt))
      break;
    momentum (grid, u, v, un, vn);
    foreach (grid) {
      u(0,0) += un(0,0)*dt;
      v(0,0) += vn(0,0)*dt;
    }
    boundary_u (grid, u, v);
    projection (grid, u, v, p, un, vn, ke);
    boundary_u (grid, u, v);
    ke_psi (grid, u, v);
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end (grid);
  cend = clock ();
  double cpu = ((double) (cend - cstart))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (N*N*(double)i/cpu));

  free_grid (grid);
}
