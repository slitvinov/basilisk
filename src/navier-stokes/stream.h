// streamfunction-vorticity Navier-Stokes solver

#include "advection.h"
#include "poisson.h"

#define omega f
scalar psi[];
mgstats mgpsi;  // statistics of the Poisson solver

// the domain boundary is a streamline i.e. psi = 0 on the boundary
psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);

event streamfunction (i++, last)
{
  // streamfunction from omega
  mgpsi = poisson (psi, omega);
  // velocity from streamfunction
  trash ({u});
  foreach_face(x)
    u.x[] = (psi[0,-1] + psi[-1,-1] - psi[0,1] - psi[-1,1])/(4.*Delta);
  foreach_face(y)
    u.y[] = (psi[1,0] + psi[1,-1] - psi[-1,0] - psi[-1,-1])/(4.*Delta);
  boundary_normal ({u});
  boundary_tangent ({u});
}
