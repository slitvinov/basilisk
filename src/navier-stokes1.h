#include "utils.h"
#include "mg.h"

scalar p[];
vector u[];
tensor S[];

// Default parameters
// Viscosity
double NU = 0.;
// Maximum number of multigrid iterations
int NITERMAX = 100;
// Tolerance on maximum divergence
double TOLERANCE = 1e-3;
// user-provided functions
void parameters (void);
void init       (void);
void end        (void);

static void boundary_u (vector u)
{
#if 0
  // slip-walls by default
  u.x[right]  = 0.;
  // u.x[left]   = 0.;
  u.y[top]    = 0.;
  // u.y[bottom] = 0.;

  boundary ((scalar *){u});
  /* fixme */
  foreach_boundary (left,true)
    u.x[] = 0.;
  foreach_boundary (bottom,true)
    u.y[] = 0.;
#else
  foreach_boundary (right,false)
    u.x[ghost] = 0.;
  foreach_boundary (top,false)
    u.y[ghost] = 0.;
  foreach_boundary (left,false)
    u.x[] = 0.;
  foreach_boundary (bottom,false)
    u.y[] = 0.;
  boundary_tangent ({u});
#endif
}

double timestep()
{
  double dtmax = DT/CFL;
  foreach_face (reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = delta/fabs(u.x[]);
      if (dt < dtmax) dtmax = dt;
    }
  return dtmax*CFL;
}

void stresses()
{
  trash ({S});
  foreach()
    foreach_dimension()
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*NU*(u.x[1,0] - u.x[])/delta;
  foreach_vertex()
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      NU*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/delta;
  boundary ({S.x.x, S.y.y});
}

void advance (double dt)
{
  foreach_face()
    u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/delta;
  boundary_u (u);
}

void relax (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - delta*delta*b[])/4.;
}

double residual (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
  foreach (reduction(max:maxres)) {
    res[] = b[] + (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/(delta*delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

void projection (vector u, scalar p, 
		 scalar div, scalar res, scalar dp)
{
  double sum = 0.;
  foreach (reduction(+:sum)) {
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/delta;
    sum += div[];
  }
  double maxres = residual (p, div, res);
  int i;
  for (i = 0; i < NITERMAX && (i < 1 || maxres > TOLERANCE); i++) {
    mg_cycle (p, res, dp,
	      relax, boundary_level,
	      4, 0);
    boundary_level ({p}, depth());
    maxres = residual (p, div, res);
  }
  if (i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  sum: %g\n", 
	     NITERMAX, sum);
  foreach_face()
    u.x[] -= (p[] - p[-1,0])/delta;
  boundary_u (u);
}

void run (void)
{
  parameters();
  init_grid (N);

  S.x.y = S.y.x; // fixme: the tensor is symmetric
  // staggering for u.x, u.y, S.x.y
  u.x.d.x = u.y.d.y = S.x.y.d.x = S.x.y.d.y = -1;
  foreach_face()
    u.x[] = 0.;
  foreach()
    p[] = 0.;
  init();
  boundary_u (u);
  boundary ({p});

  projection (u, p, S.x.y, S.y.y, S.x.x);

  timer start = timer_start();
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep());
    stresses();
    advance (dt);
    projection (u, p, S.x.y, S.y.y, S.x.x);
    i++; t = tnext;
  }
  end();
  timer_print (start, i, -1);

  free_grid();
}
