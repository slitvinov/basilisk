// streamfunction-vorticity Navier-Stokes solver

#include "utils.h"
#include "timestep.h"
#include "bcg.h"
#include "poisson.h"

scalar omega[], psi[];
vector u[];

// Default parameters
// gradient
double (* gradient) (double, double, double) = NULL; // centered
// user-provided functions
void parameters     (void);
void init           (void);

double dt = 0.; // timestep
mgstats mgpsi;  // statistics of the Poisson solver

// no flow through boundaries
u.x[right]  = 0.;
u.x[left]   = 0.;
u.y[top]    = 0.;
u.y[bottom] = 0.;

// the domain boundary is a streamline i.e. psi = 0 on the boundary
psi[right]  = - psi[];
psi[left]   = - psi[];
psi[top]    = - psi[];
psi[bottom] = - psi[];

void run (void)
{
  CFL = 0.8; // otherwise it's 0.5

  parameters();
  init_grid (N);

  u.x.d.x = u.y.d.y = -1; // staggering for u.x, u.y
  omega.gradient = gradient;
  foreach_face()
    u.x[] = 0.;
  foreach()
    omega[] = psi[] = 0.;
  init();
  boundary ({omega, psi});
  boundary_mac ({u});

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    // streamfunction from omega
    mgpsi = poisson (psi, omega);
    // velocity from streamfunction
    foreach_face(x)
      u.x[] = (psi[0,-1] + psi[-1,-1] - psi[0,1] - psi[-1,1])/(4.*Delta);
    foreach_face(y)
      u.y[] = (psi[1,0] + psi[1,-1] - psi[-1,0] - psi[-1,-1])/(4.*Delta);
    boundary_mac ({u});
    // advection of vorticity
    dt = dtnext (t, timestep (u));
    vector flux[];
    fluxes_upwind_bcg (omega, u, flux, dt);
    foreach(reduction(+:tnc)) {
      omega[] += dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/Delta;
      tnc++;
    }
    boundary ({omega});
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
