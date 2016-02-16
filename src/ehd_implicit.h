/**
# Electrohydrodynamics ohmic conduction

This function approximates implicitly the EHD equation set given by
the electric Poisson equation and the ohmic conduction term of the
charge density conservation (the advection term is computed elsewhere
using a tracer advection scheme),
$$
\nabla \cdot( \epsilon \nabla  \phi^{n+1}) = -\rho_e^{n+1} \quad \mbox{and}
\quad (\rho_e^{n+1}-\rho_e^n) = \Delta t\nabla \cdot (K\nabla \phi^{n+1})
$$ 
where $\rho_e$ is the charge density, and $\phi$ is the electric potential.
$K$ and $\epsilon$ are the conductivity and permittivity, respectively.

Substituting the Poisson equation into the conservation equation gives
the following time-implicit scheme,
$$
\nabla \cdot [(K \, \Delta t + \epsilon) \nabla \phi^{n+1}] = -\rho_e^{n}
$$
We need the Poisson solver and the timestep *dt*. */

#include "poisson.h"
extern double dt;

/**
The electric potential and the volume charge density are scalars while
the permittivity and conductivity are face vectors. In *mgphi* we will
store the statistics for the multigrid resolution of the electric
poisson equation. */

scalar phi[], rhoe[];
face vector epsilon[], K[];
mgstats mgphi;

/**
In case the variables were not initialized, the event default set the
initial values. */

event defaults (i = 0) {

  /**
  The coarsen/refine attributes of the charge density are those of a tracer
  otherwise the conservation is not guaranteed. */

#if QUADTREE
  rhoe.coarsen = coarsen_volume_average;
  rhoe.refine = refine_linear;
#endif

  /**
  By default, both the conductivity and permittivity are unity and the
  charge density and electric fields are zero. */
  
  foreach_face()
    epsilon.x[] = K.x[] = fm.x[];
  foreach()
    phi[] = rhoe[] = 0.;
  boundary ({phi,rhoe});
}

event ohmic_conduction (i++) {
  scalar rhs[];
  
  /**
  The r.h.s of the poisson equation has to to be set to $-\rho_e^n$, */ 
  
  foreach()
    rhs[] = - rhoe[]*cm[];
  if (K.x.i) {
    
    /**
    The face vector $(K dt +\epsilon)$ is constructed first, */

    face vector f[];
    foreach_face()
      f.x[] = K.x[]*dt + epsilon.x[];
    boundary_flux ({f});

    /**
    The poisson equation is solved to obtain $\phi^{n+1}$, */

    mgphi = poisson (phi, rhs, f);

    /**
    Finally, $\rho_e^{n+1}$ is calculated by computing 
    $\rho_e^{n+1} = -\nabla \cdot (\epsilon \nabla \phi^{n+1})$. */

#if QUADTREE
    foreach_face()
      f.x[] = epsilon.x[]*(phi[] - phi[-1])/Delta;
    boundary_flux ({f});
    foreach() {
      rhoe[] = 0.;
      foreach_dimension()
	rhoe[] -= f.x[1] - f.x[];
      rhoe[] /= cm[]*Delta;
    }
#else // Cartesian
    foreach() {
      rhoe[] = 0.;
      foreach_dimension()
	rhoe[] -= (epsilon.x[1]*(phi[1] - phi[]) -
		   epsilon.x[]*(phi[] - phi[-1]));
      rhoe[] /= cm[]*sq(Delta);
    }
#endif
    boundary ({rhoe});
  }

  /**
  In the absence of conductivity, the electric potential only varies 
  if the electrical permittivity varies, */
  
  else
    mgphi = poisson (phi, rhs, epsilon);
}
