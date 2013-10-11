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
void end        (void);
// timestep
double dt = 0;
mgstats mgp;  // statistics of the Poisson solver

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
      S.x.x[] = - sq(u.x[] + u.x[1,0])/4. + 2.*NU*(u.x[1,0] - u.x[])/Delta;
  foreach_vertex()
    S.x.y[] = 
      - (u.x[] + u.x[0,-1])*(u.y[] + u.y[-1,0])/4. +
      NU*(u.x[] - u.x[0,-1] + u.y[] - u.y[-1,0])/Delta;
  boundary ({S.x.x, S.y.y});

  // update
  foreach_face()
    u.x[] += dt*(S.x.x[] - S.x.x[-1,0] + S.x.y[0,1] - S.x.y[])/Delta;
  boundary_mac ({u});
}

void projection (vector u, scalar p)
{
  scalar div[];
  foreach()
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/Delta;
  mgp = poisson (p, div);
  foreach_face()
    u.x[] -= (p[] - p[-1,0])/Delta;
  boundary_mac ({u});
}

event defaults (i = 0)
{
  // staggering for u.x, u.y
  u.x.d.x = u.y.d.y = -1;
  foreach_face()
    u.x[] = 0.;
  boundary_mac ({u});
  foreach()
    p[] = 0.;
  boundary ({p});
}

event init (i = 0)
{
  boundary_mac ({u});
  boundary ({p});
  projection (u, p);
}

void run (void)
{
  parameters();
  init_grid (N);

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
