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

// slip-walls by default
u.x[right]  = 0.;
// u.x[left]   = 0.;
u.y[top]    = 0.;
// u.y[bottom] = 0.;

static void boundary_uv (vector u)
{
  boundary ((scalar *){u});
  /* fixme */
  foreach_boundary (left,true)
    u.x[] = 0.;
  foreach_boundary (bottom,true)
    u.y[] = 0.;
}

double timestep ()
{
  double dtmax = DT/CFL;
  foreach(reduction(min:dtmax)) {
    double dx = L0*delta;
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = dx/fabs(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  }
  return dtmax*CFL;
}

#define sq(x) ((x)*(x))

void stresses()
{
  S.x.y = S.y.x; // fixme: the tensor is symmetric
  foreach() {
    foreach_dimension()
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*NU*(u.x[1,0] - u.x[])/delta;
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      NU*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/delta;
  }
  foreach_boundary (left, false)
    S.x.x[ghost] = - sq(u.x[-1,0] + u.x[])/4. + 2.*NU*(u.x[] - u.x[-1,0])/delta;
  foreach_boundary (top, false)
    S.x.y[ghost] = - (u.x[0,1] + u.x[])*(u.y[0,1] + u.y[-1,1])/4. +
      NU*(u.x[0,1] - u.x[] + u.y[0,1] - u.y[-1,1])/delta;
  foreach_boundary (right, false)
    S.x.y[ghost] = - (u.x[1,0] + u.x[1,-1])*(u.y[1,0] + u.y[])/4. +
      NU*(u.x[1,0] - u.x[1,-1] + u.y[1,0] - u.y[])/delta;
  foreach_boundary (bottom, false)
    S.y.y[ghost] = - sq(u.y[0,-1] + u.y[])/4. + 2.*NU*(u.y[] - u.y[0,-1])/delta;
}

void advance (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/delta;
}

void relax (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - L0*L0*delta*delta*b[])/4.;
}

double residual (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
  foreach(reduction(max:maxres)) {
    res[] = b[] + (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])
      /(L0*L0*delta*delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

void projection (vector u, scalar p, 
		 scalar div, scalar res, scalar dp)
{
  double sum = 0.;
  foreach(reduction(+:sum)) {
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/(L0*delta);
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
}

void run (void)
{
  parameters ();
  init_grid (N);

  foreach()
    u.x[] = u.y[] = p[] = 0.;
  init ();
  boundary_uv (u);
  boundary ({p});

  projection (u, p, S.x.y, S.y.y, S.x.x);
  boundary_uv (u);

  timer start = timer_start();
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep ());
    stresses ();
    advance (dt);
    boundary_uv (u);
    projection (u, p, S.x.y, S.y.y, S.x.x);
    boundary_uv (u);
    i++; t = tnext;
  }
  end ();
  timer_print (start, i, -1);

  free_grid();
}
