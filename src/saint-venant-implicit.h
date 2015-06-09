/**
# A semi-implicit Saint-Venant solver

We solve the [Saint-Venant equations](saint-venant.h) semi-implicitly
to lift the timestep restriction due to gravity waves. This is just a
particular application of the "all Mach" semi-implicit solver. We will
use the Bell-Collela-Glaz advection scheme to transport mass and
momentum. */

#include "all-mach.h"
#include "tracer.h"

/**
In addition to the momentum $\mathbf{q}=h\mathbf{u}$ defined by the
all Mach solver, we will need the fluid depth *h* (i.e. the density
$\rho$) and the topography *zb*. Both the momentum $\mathbf{q}$ and
mass $h$ are advected tracers. The acceleration of gravity is
*G*. *dry* is the fluid depth below which a cell is considered
"dry". */

scalar h[], zb[], * tracers = {q,h};
double G = 1.;
double dry = 1e-10;

/**
We need fields to store the (varying) state field $\rho c^2$ and the
acceleration $\mathbf{a}$. */

scalar rhoc2v[];
face vector av[];

event defaults (i = 0) {
  rho = h;
  a = av;
  rhoc2 = rhoc2v;

  /**
  We set limiting for *q*, *h* and *zb*. */
  
  for (scalar s in {q,h,zb})
    s.gradient = minmod2;

  /**
  As well as default values. */

  foreach() {
    h[] = 1.;
    zb[] = 0.;
  }
  boundary ({h,zb});
  
#if QUADTREE
  for (scalar s in {q,h,zb}) {
    s.refine = s.prolongation = refine_linear;
    s.coarsen = coarsen_volume_average;
  }
#endif
}

/**
We apply boundary conditions after user initialisation. The boundary
conditions for $\rho=h$ and $\mathbf{q}$ are already applied by the
[all Mach solver](all-mach.h). */

event init (i = 0) {
  boundary ({zb});
}

/**
The properties of the fluid are given by the "Saint-Venant equation of
state" i.e. $p = gh^2/2$, $c^2 = gh$. */

event properties (i++) {

  /**
  We first make sure that the depth is (slightly) larger than zero. */
  
  foreach()
    if (h[] < dry)
      h[] = dry;
  boundary ({h});

  foreach() {
    rhoc2v[] = G*sq(h[]);
    ps[] = rhoc2v[]/2.;
  }
}

/**
The acceleration due to the topography is $- g\nabla z_b$. */

event acceleration (i++) {
  foreach_face()
    av.x[] += - G*(zb[] - zb[-1])/Delta;
}
