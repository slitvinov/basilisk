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
The scheme is stable for a CFL smaller than 0.5.

On quadtrees, we need to setup the appropriate prolongation and
refinement functions for the volume fraction fields. */

event defaults (i = 0)
{
  CFL = 0.5;
#if QUADTREE
  for (scalar c in interfaces) {
    c.prolongation = fraction_prolongation;
    c.refine = fraction_refine;
  }
#endif
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

/**
We first reconstruct the interface normal $\mathbf{n}$ and the
intercept $\alpha$ for each cell. Then we go through each (vertical)
face of the grid. */

  reconstruction (c, n, alpha);
  foreach_face(x) {

/**
To compute the volume fraction flux, we check the sign of the velocity
component normal to the face and compute the index `i` of the
corresponding *upwind* cell (either 0 or -1). */

    double un = uf.x[]*dt/Delta, s = sign(un);
    int i = -(s + 1.)/2.;

/**
If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
volume fraction flux through the face of the cell is given by the dark
area in the figure below. The corresponding volume fraction can be
computed using the `rectangle_fraction()` function.

![Volume fraction flux](figures/flux.png)

When the upwind cell is entirely full or empty we can avoid this
computation. */

    double cf = (c[i,0] <= 0. || c[i,0] >= 1.) ? c[i,0] :
      rectangle_fraction (- s*n.x[i,0], n.y[i,0], alpha[i,0],
			  -0.5, -0.5, s*un - 0.5, 0.5);
/**
Once we have the upwind volume fraction, the volume fraction flux
through the face is simply: */

    flux[] = uf.x[]*cf;
  }

/**
On quadtree grids, we need to make sure that the fluxes match at
fine/coarse cell boundaries i.e. we need to *restrict* the fluxes from
fine cells to coarse cells. This is what is usually done, for all
dimensions, by the `boundary_normal()` function. Here, we only need to
do it for a single dimension (x). */

#if QUADTREE
  foreach_halo_fine_to_coarse() {
    flux[] = (fine(flux,0,0) + fine(flux,0,1))/2.;
    if (is_leaf (neighbor(1,0)))
      flux[1,0] = (fine(flux,2,0) + fine(flux,2,1))/2.;
  }
#endif

/**
Once we have computed the fluxes on all faces, we can update the
volume fraction field according to the one-dimensional advection
equation
$$
\partial_tc = \mathbf{u}_f\cdot\nabla_x c + c\nabla_x\cdot\mathbf{u}_f
$$
The first term is computed using the fluxes. The second term -- which is
non-zero for the one-dimensional velocity field -- is approximated using
a centered volume fraction field `cc` which will be defined below. */

  foreach()
    c[] += dt*(flux[] - flux[1,0] + cc[]*(uf.x[1,0] - uf.x[]))/Delta;
  boundary ({c});
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

event vof (i = 1; i++)
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

    void (* sweep[2]) (scalar, scalar) = {sweep_x, sweep_y};
    for (int d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc);
  }
}
