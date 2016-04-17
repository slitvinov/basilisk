/**
# Volume-Of-Fluid advection

We want to approximate the solution of the advection equations
$$
\partial_tc_i + \mathbf{u}_f\cdot\nabla c_i = 0
$$
where $c_i$ are volume fraction fields describing sharp interfaces.

This can be done using a conservative, non-diffusive geometric VOF
scheme.

We will need basic functions for volume fraction computations. */

#include "fractions.h"

/**
The list of volume fraction fields `interfaces`, will be provided by
the user.

The face velocity field `uf` will be defined by a solver as well
as the timestep. */

extern scalar * interfaces;
extern face vector uf;
extern double dt;

/**
On trees, we need to setup the appropriate prolongation and
refinement functions for the volume fraction fields. */

event defaults (i = 0)
{
#if TREE
  for (scalar c in interfaces)
    c.refine = c.prolongation = fraction_refine;
#endif
}

/**
We need to make sure that the CFL is smaller than 0.5 to ensure
stability of the VOF scheme. */

event stability (i = 0) {
  if (CFL > 0.5)
    CFL = 0.5;
}

/**
## One-dimensional advection

The simplest way to implement a multi-dimensional VOF advection scheme
is to use dimension-splitting i.e. advect the field along each
dimension successively using a one-dimensional scheme.

We implement the one-dimensional scheme along the x-dimension and use
the [foreach_dimension()](/Basilisk C#foreach_dimension) operator to
automatically derive the corresponding functions along the other
dimensions. */

foreach_dimension()
static void sweep_x (scalar c, scalar cc)
{
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;

  /**
  We first reconstruct the interface normal $\mathbf{n}$ and the
  intercept $\alpha$ for each cell. Then we go through each (vertical)
  face of the grid. */

  reconstruction (c, n, alpha);

  foreach_face(x, reduction (max:cfl)) {

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

    double un = uf.x[]*dt/(Delta*fm.x[]), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

    if (un*fm.x[]*s/cm[] > cfl)
      cfl = un*fm.x[]*s/cm[];

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.
    
    ![Volume fraction flux](figures/flux.svg)
    
    When the upwind cell is entirely full or empty we can avoid this
    computation. */

    double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5});
    
    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    flux[] = cf*uf.x[];
  }

  /**
  On tree grids, we need to make sure that the fluxes match at
  fine/coarse cell boundaries i.e. we need to *restrict* the fluxes from
  fine cells to coarse cells. This is what is usually done, for all
  dimensions, by the `boundary_flux()` function. Here, we only need to
  do it for a single dimension (x). */

#if TREE
  for (int l = depth() - 1; l >= 0; l--)
    foreach_halo (prolongation, l) {
#if dimension == 1
      if (is_refined (neighbor(-1)))
	flux[] = fine(flux,0);
      if (is_refined (neighbor(1)))
	flux[1] = fine(flux,2);
#elif dimension == 2
      if (is_refined (neighbor(-1)))
	flux[] = (fine(flux,0,0) + fine(flux,0,1))/2.;
      if (is_refined (neighbor(1)))
	flux[1] = (fine(flux,2,0) + fine(flux,2,1))/2.;
#else // dimension == 3
      if (is_refined (neighbor(-1)))
	flux[] = (fine(flux,0,0,0) + fine(flux,0,1,0) +
		  fine(flux,0,0,1) + fine(flux,0,1,1))/4.;
      if (is_refined (neighbor(1)))
	flux[1] = (fine(flux,2,0,0) + fine(flux,2,1,0) +
		   fine(flux,2,0,1) + fine(flux,2,1,1))/4.;
#endif
    }
#endif

  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  /**
  Once we have computed the fluxes on all faces, we can update the
  volume fraction field according to the one-dimensional advection
  equation
  $$
  \partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  $$
  The first term is computed using the fluxes. The second term -- which is
  non-zero for the one-dimensional velocity field -- is approximated using
  a centered volume fraction field `cc` which will be defined below. */

  foreach()
    c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  boundary ({c});
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

event vof (i++)
{
  for (scalar c in interfaces) {

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[];
    foreach()
      cc[] = (c[] > 0.5);

    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    boundary ({c});
    for (d = 0; d < dimension; d++)
      sweep[(i + d) % dimension] (c, cc);
  }
}
