#include "utils.h"
#include "timestep.h"
#include "poisson.h"

scalar p[];
vector u[];

// Default parameters
// Viscosity
double NU = 0.;
// user-provided functions
void parameters (void);
void init       (void);
void end        (void);
// timestep
double dt = 0;

// default is no flow through boundaries
u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

void advance (double dt)
{
  // stresses
  symmetric tensor S[];
  // staggering for S.x.y
  S.x.y.d.x = S.x.y.d.y = -1;

  trash ({S});
  foreach()
    foreach_dimension()
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*NU*(u.x[1,0] - u.x[])/delta;
  foreach_vertex()
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      NU*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/delta;
  boundary ({S.x.x, S.y.y});

  // update
  foreach_face()
    u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/delta;
  boundary_mac ({u});
}

void projection (vector u, scalar p)
{
  scalar div[], res[], dp[];
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
    // fixme: use a generic Poisson solver instead
    // need to find a solution for homogeneous BCs
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
  boundary_mac ({u});
}

void run (void)
{
  parameters();
  init_grid (N);

  // staggering for u.x, u.y
  u.x.d.x = u.y.d.y = -1;
  foreach_face()
    u.x[] = 0.;
  foreach()
    p[] = 0.;
  init();
  boundary_mac ({u});
  boundary ({p});

  projection (u, p);

  timer start = timer_start();
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    double dt = dtnext (t, timestep (u));
    advance (dt);
    projection (u, p);
    i++; t = tnext;
  }
  end();
  timer_print (start, i, -1);

  free_grid();
}
