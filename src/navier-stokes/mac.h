#include "utils.h"
#include "timestep.h"
#include "poisson.h"

scalar p[];
staggered vector u[];

// Default parameters
// Viscosity
double NU = 0.;
// specific volume (default is unity)
staggered vector alpha;
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

event defaults (i = 0)
{
  foreach_face()
    u.x[] = 0.;
  foreach()
    p[] = 0.;
  boundary ({p,u});
}

event init (i = 0)
{
  boundary ({p,u});
}

event advance (i++)
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
}

event projection (i++)
{
  scalar div[];
  boundary ((scalar *){u});
  foreach()
    div[] = (u.x[1,0] - u.x[] + u.y[0,1] - u.y[])/Delta;
  mgp = poisson (p, div, alpha);
  if (alpha.x)
    foreach_face()
      u.x[] -= alpha.x[]*(p[] - p[-1,0])/Delta;
  else
    foreach_face()
      u.x[] -= (p[] - p[-1,0])/Delta;
  boundary ((scalar *){u});

  dt = dtnext (t, timestep (u));
}

void run (void)
{
  parameters();
  init_grid (N);

  timer start = timer_start();
  double t = 0;
  int i = 0;
  while (events (i, t)) {
    i++; t = tnext;
  }
  end();
  timer_print (start, i, -1);

  free_grid();
}
