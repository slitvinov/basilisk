/**
# Bell-Collela-Glaz advection scheme

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89). Given a centered
scalar field *f*, a face vector field *u*, a timestep *dt* and a
source term field *src*, it fills the face vector field *flux* with
the components of the advection fluxes of *f*. */

void tracer_fluxes (const scalar f, 
		    const face vector u,
		    face vector flux,
		    double dt,
		    (const) scalar src)
{

  /**
  We first compute the cell-centered gradient of *f* in a locally-allocated
  vector field. */
  
  vector g[];
  gradients ({f}, {g});

  /**
  For each face, the flux is composed of two parts... */

  foreach_face() {

    /**
    A normal component... */

    double un = dt*u.x[]/Delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i,0] + src[i,0]*dt/2. +
      s*min(1., 1. - s*un)*g.x[i,0]*Delta/2.;

    /**
    and a tangential component... */

    double vn = u.y[i,0] + u.y[i,1];
    double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
    f2 -= dt*vn*fyy/(4.*Delta);
    flux.x[] = f2*u.x[];
  }

  /**
  Boundary conditions ensure the consistency of fluxes across
  variable-resolution boundaries (on adaptive meshes). */

  boundary_flux ({flux});
}

/**
The function below uses the *tracer_fluxes* function to integrate the
advection equation, using an explicit scheme with timestep *dt*, for
each tracer in the list. */

struct Advection {
  scalar * tracers;
  face vector u;
  double dt;
  scalar * src; // optional
};

void advection (struct Advection p)
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  for (f,src in p.tracers,lsrc) {
    vector flux[];
    tracer_fluxes (f, p.u, flux, p.dt, src);
    foreach()
      f[] += p.dt*(flux.x[] - flux.x[1,0] + flux.y[] - flux.y[0,1])/Delta;
    boundary ({f});
  }

  if (!p.src)
    free (lsrc);
}
