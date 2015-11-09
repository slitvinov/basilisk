/**
# Surface tension

We will need to compute the curvature of the interface, using its
Volume-Of-Fluid description. */

#include "curvature.h"

/**
The surface tension $\sigma$ and interface curvature $\kappa$ will be
associated to each VOF tracer. This is done easily by adding the
following [field attributes](/Basilisk C#field-attributes). */

attribute {
  double sigma;
  scalar kappa;
}

/**
Surface tension is a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {
  if (is_constant(a.x))
    a = new face vector;
}

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  We then consider each VOF interface with an associated value of
  $\sigma$ different from zero and set the maximum timestep. */
  
  for (scalar c in interfaces)
    if (c.sigma) {
      double dt = sqrt (rhom*cube(dmin)/(pi*c.sigma));
      if (dt < dtmax)
	dtmax = dt;
    }
}

/**
## Surface tension term

The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  
  /**
  We check for all VOF interfaces for which $\sigma$ is non-zero. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;
  for (scalar c in interfaces)
    if (c.sigma) {
      list = list_add (list, c);

      /**
      To avoid undeterminations due to round-off errors, we remove
      values of the volume fraction larger than one or smaller than
      zero. */

      foreach()
	c[] = clamp (c[], 0, 1);
      boundary ({c});
      
      /**
      We also allocate the curvature fields (if this has not been done
      already) and update the values using the height-function
      curvature calculation. */

      if (!c.kappa.i)
	c.kappa = new scalar;
      curvature (c, c.kappa);
    }

  /**
  Finally, for each interface for which $\sigma$ is non-zero, we
  compute the surface tension acceleration
  $$
  \sigma\kappa\mathbf{n}\delta_s/\rho \approx \alpha\sigma\kappa\nabla c
  $$ 
  */

  face vector st = a;
  foreach_face()
    for (scalar c in list)
      if (c[] != c[-1]) {
	scalar kappa = c.kappa;
	
	/**
	We need to compute the curvature *kf* on the face, using its
	values at the center of the cell. If both cells are cut by the
	interface, we take the average of both curvatures, otherwise
	we take the single value in the cell cut by the interface
	(*nodata* is a very large positive value). */

	double kf = 
	  (kappa[] != nodata && kappa[-1] != nodata) ? 
	  (kappa[] + kappa[-1])/2. : 
	  min(kappa[], kappa[-1]);

	st.x[] += alpha.x[]/fm.x[]*c.sigma*kf*(c[] - c[-1])/Delta;
      }

  /**
  Finally we free the list of interfacial volume fractions. */

  free (list);
}
