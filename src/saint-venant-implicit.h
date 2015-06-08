/**
# A semi-implicit Saint-Venant solver

We solve the [Saint-Venant equations](saint-venant.h) semi-implicitly
to lift the timestep restriction due to gravity waves. This is just a
particular application of the "all Mach" semi-implicit solver. */

#include "all-mach.h"

/**
In addition to the momentum $\mathbf{q}=h\mathbf{u}$ defined by the
all Mach solver, we will need the fluid depth *h* (i.e. the density
$\rho$) and the topography *zb*. The acceleration of gravity is
*G*. *dry* is the fluid depth below which a cell is considered
"dry". */

scalar h[], zb[];
double G = 1.;
double dry = 1e-10;

/**
We need fields to store the (varying) state fields $\alpha=1/\rho$,
$\rho c^2$ and the acceleration $\mathbf{a}$. */

face vector alphav[], av[];
scalar rhoc2v[];

event defaults (i = 0) {
  alpha = alphav;
  a = av;
  rhoc2 = rhoc2v;

  /**
  We set limiting for *h* and *zb*. */
  
  for (scalar s in {h,zb})
    s.gradient = minmod2;

  /**
  As well as default values. */

  foreach() {
    h[] = 1.;
    zb[] = 0.;
  }
  boundary ({h,zb});
  
#if QUADTREE
  for (scalar s in {h,zb}) {
    s.refine = s.prolongation = refine_linear;
    s.coarsen = coarsen_volume_average;
  }
#endif
}

/**
We make sure that *h* is larger than *dry* after user initialisation. */

event init (i = 0) {
  foreach()
    if (h[] < dry)
      h[] = dry;
  boundary ({h});
}

/**
The water depth is advected by the face centered velocity field. */

event tracer_advection (i++) {
  advection ({h}, uf, dt);
  foreach()
    if (h[] < dry)
      h[] = dry;
  boundary ({h});
}

/**
The properties of the fluid are given by the "Saint-Venant equation of
state" i.e. $p = gh^2/$, $c^2 = gh$, $\alpha = 1/h$. */

event properties (i++) {
  foreach() {
    rhoc2v[] = G*sq(h[]);
    ps[] = rhoc2v[]/2.;
  }
  foreach_face()
    alphav.x[] = 2./(h[] + h[-1]);
}

/**
The acceleration is $- g\nabla z_b$. */

event acceleration (i++) {
  foreach_face()
    av.x[] += - G*(zb[] - zb[-1])/Delta;
}
