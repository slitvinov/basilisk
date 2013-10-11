#include "utils.h"

vector u[], un[];
scalar h[], hn[], b[];
scalar ke[], psi[];

// Default parameters
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;
// Viscosity
double NU = 0.;
// user-provided functions
void parameters (void);

void boundary_uv (vector u)
{
  /* slip walls (symmetry) by default */
  foreach_boundary (right, true)
    u.x[ghost] = 0.;
  foreach_boundary (right, true)
    u.y[ghost] = u.y[];
  foreach_boundary (left, true) {
    u.x[ghost] = - u.x[1,0];
    u.x[] = 0.;
  }
  foreach_boundary (left, true)
    u.y[ghost] = u.y[];
  foreach_boundary (top, true)
    u.x[ghost] = u.x[];
  foreach_boundary (top, true)
    u.y[ghost] = 0.;
  foreach_boundary (bottom, true)
    u.x[ghost] = u.x[];
  foreach_boundary (bottom, true) {
    u.y[ghost] = - u.y[0,1];
    u.y[] = 0.;
  }
}

void advection_centered (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((f[] + f[-1,0])*u.x[] - 
	    (f[] + f[1,0])*u.x[1,0] +
	    (f[] + f[0,-1])*u.y[] - 
	    (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
}

void advection_upwind (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] - 
	    (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
	    (u.y[] < 0. ? f[] : f[0,-1])*u.y[] - 
	    (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
}
    
double timestep (void)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach(reduction(min:dtmax)) {
    Delta *= Delta;
    if (h[] > 0.) {
      double dt = Delta/(G*h[]);
      if (dt < dtmax) dtmax = dt;
    }
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = Delta/sq(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  }
  return sqrt (dtmax)*CFL;
}

void momentum (vector u, scalar h, vector du)
{
  struct { double x, y; } f = {1.,-1.};
  foreach() {
    double g = G*(h[] + b[]) + ke[];
    foreach_dimension() {
      double psiu = (psi[] + psi[0,1])/2.;
      du.x[] = 
	- (g - G*(h[-1,0] + b[-1,0]) - ke[-1,0])/Delta
	+ f.x*(psiu + F0)*(u.y[] + u.y[0,1] + u.y[-1,0] + u.y[-1,1])/4.
	+ NU*(u.x[1,0] + u.x[0,1] + u.x[-1,0] + u.x[0,-1] - 4.*u.x[])/sq(Delta);
    }
  }
}

void ke_psi (vector u)
{
  foreach() {
#if 1
    ke[] = (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.;
#else
    double uc = u.x[]*u.x[] + u.x[1,0]*u.x[1,0];
    double vc = u.y[]*u.y[] + u.y[0,1]*u.y[0,1];
    ke[] = (uc + vc)/4.;
#endif
    psi[] = (u.y[] - u.y[-1,0] + u.x[0,-1] - u.x[])/Delta;
  }
  foreach_boundary (right, false)
    psi[1,0] = (u.y[1,0] - u.y[] + u.x[1,-1] - u.x[1,0])/Delta;
  foreach_boundary (left, false)
    ke[-1,0] = (sq(u.x[-1,0] + u.x[]) + sq(u.y[-1,0] + u.y[-1,1]))/8.;
  foreach_boundary (top, false)
    psi[0,1] = (u.y[0,1] - u.y[-1,1] + u.x[] - u.x[0,1])/Delta;
  foreach_boundary (bottom, false)
    ke[0,-1] = (sq(u.x[0,-1] + u.x[1,-1]) + sq(u.y[0,-1] + u.y[]))/8.;
}

void advance (double t, scalar * f, scalar * df)
{
  vector u = {f[0], f[1]}, du = {df[0], df[1]};
  scalar h = f[2], dh = df[2];

  advection_centered (h, u, dh);
  momentum (u, h, du);
}

void update (double t, scalar * f)
{
  vector u = {f[0], f[1]};
  scalar h = f[2];
  boundary ({h});
  boundary_uv (u);
  ke_psi (u);
}

event defaults (i = 0)
{
  foreach() {
    b[] = u.x[] = u.y[] = 0.;
    h[] = 1.;
  }
  boundary ({b,h});
  boundary_uv (u);
  ke_psi (u);
}

event init (i = 0)
{
  boundary ({b,h});
  boundary_uv (u);
  ke_psi (u);
}

void run (void)
{
  parameters ();
  init_grid (N);

  timer start = timer_start();
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep ());
#if 1
    advection_centered (h, u, hn);
    foreach()
      h[] += hn[]*dt;
    boundary ({h});
    momentum (u, h, un);
    foreach()
      foreach_dimension()
        u.x[] += un.x[]*dt;
    boundary_uv (u);
    ke_psi (u);
#else /* unstable! */
    scalar f[3] = { u, v, h };
    scalar df[2][3] = {{ un,  vn,  hn },
		       { un1, vn1, hn1 }};
    runge_kutta (2, t, dt, 3, f, df, advance, update);
#endif
    i++; t = tnext;
  }
  timer_print (start, i, -1);

  free_grid ();
}
