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

event defaults (i = 0) {
  
  /**
  Surface tension is a source term in the right-hand-side of the
  evolution equation for the velocity of the [centered Navier--Stokes
  solver](navier-stokes/centered.h) i.e. it is an acceleration. If
  necessary, we allocate a new vector field to store it. */

  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }

  /**
  Each interface for which $\sigma$ is not zero needs a new field to
  store the curvature. */
  
  for (scalar c in interfaces)
    if (c.sigma && !c.kappa.i) {
      scalar kappa = new_scalar ("kappa");
      foreach()
	kappa[] = 0.;
      boundary ({kappa});
      c.kappa = kappa;
    }
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
      We update the values of $\kappa$ using the height-function
      curvature calculation. */

      assert (c.kappa.i);
      curvature (c, c.kappa);
    }

  /**
  On trees we need to make sure that the volume fraction gradient
  is computed exactly like the pressure gradient. This is necessary to
  ensure well-balancing of the pressure gradient and surface tension
  term. To do so, we apply the same prolongation to the volume
  fraction field as applied to the pressure field. */
  
#if TREE
  for (scalar c in list)
    c.prolongation = p.prolongation;
  boundary (list);
#endif

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
	values at the center of the cell. If both curvatures are
	defined, we take the average, otherwise we take a single
	value. If all fails we set the curvature to zero: this should
	happen only because of very pathological cases e.g. weird
	boundary conditions for the volume fraction. */

	double kf = 
	  (kappa[] < nodata && kappa[-1] < nodata) ?
	     (kappa[] + kappa[-1])/2. :
	  kappa[] < nodata ? kappa[] :
	  kappa[-1] < nodata ? kappa[-1] :
	  0.;
	
	st.x[] += alpha.x[]/fm.x[]*c.sigma*kf*(c[] - c[-1])/Delta;
      }

  /**
  On trees, we need to restore the prolongation values for the
  volume fraction field. */
  
#if TREE
  for (scalar c in list)
    c.prolongation = fraction_refine;
  boundary (list);
#endif
  
  /**
  Finally we free the list of interfacial volume fractions. */

  free (list);
}
