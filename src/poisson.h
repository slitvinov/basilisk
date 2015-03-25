/**
# Multigrid Poisson--Helmholtz solvers

We want to solve Poisson--Helmholtz equations of the general form
$$
L(a) = \nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
This can be done efficiently using a multigrid solver. 

An important aspect of Poisson--Helmholtz equations is that the
operator $L()$ is linear. This property can be used to build better
estimates of a solution by successive *corrections* to an initial
guess. If we define an approximate solution $\tilde{a}$ as
$$
\tilde{a} + da = a
$$
where $a$ is the exact (unknown) solution, using the linearity of the
operator we find that $da$ verifies
$$
L(da) = b - L(\tilde{a})
$$
where the right-hand-side is often called the *residual* of the
approximate solution $\tilde{a}$.

## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
*a*, a residual *res*, a correction field *da* and a relaxation
function *relax*, we will provide an improved guess at the end of the
cycle. */

void mg_cycle (scalar * a, scalar * res, scalar * da,
	       void (* relax) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  for (int l = minlevel; l <= depth(); l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
	for (scalar s in da)
	  s[] = (9.*coarse(s,0,0) + 
		 3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;

    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      s[] += ds[];
  }
  boundary (a);
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */

int NITERMAX = 100;
double TOLERANCE = 1e-3;

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
} mgstats;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. */

mgstats mg_solve (scalar * a, scalar * b,
		  double (* residual) (scalar * a, scalar * b, scalar * res,
				       void * data),
		  void (* relax) (scalar * da, scalar * res, int depth, 
				  void * data),
		  void * data)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (a), * res = list_clone (a);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0, 0., 0.};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in b)
      sum += s[];
  s.sum = sum;

  /**
  Here we compute the initial residual field and its maximum. */

  s.resb = s.resa = residual (a, b, res, data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.resa > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, data, 4, 0);
    s.resa = residual (a, b, res, data);
  }

  /**
  If we have reached the maximum number of iterations, we warn the user. */

  if (s.i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  res: %g sum: %g\n", 
	     NITERMAX, s.resa, s.sum);

  /**
  We deallocate the residual and correction fields and free the lists. */

  delete (res); free (res);
  delete (da);  free (da);
  
  return s;
}

/**
## Application to the Poisson--Helmholtz equation

We now apply the generic multigrid solver to the Poisson--Helmholtz equation
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
We first setup the data structure required to pass the extra
parameters $\alpha$ and $\lambda$. We define $\alpha$ as a face
vector field because we need values at the face locations
corresponding to the face gradients of field $a$. 

*alpha* and *lambda* are declared as *(const)* to indicate that the
function works also when *alpha* and *lambda* are constant vector
(resp. scalar) fields. If *tolerance* is set, it supersedes the
default *TOLERANCE* of the multigrid solver. */

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l)
    a[] = (alpha.x[1,0]*a[1,0] + alpha.x[]*a[-1,0] + 
	   alpha.y[0,1]*a[0,1] + alpha.y[]*a[0,-1] 
	   - sq(Delta)*b[])
    /(alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[] - lambda[]*sq(Delta));

#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the quadtree mesh
requires more careful consideration... */

static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1,0])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta
      - lambda[]*a[];
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + 
      ((alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[])*a[]
       - alpha.x[1,0]*a[1,0] - alpha.x[]*a[-1,0] 
       - alpha.y[0,1]*a[0,1] - alpha.y[]*a[0,-1])/sq(Delta)
      - lambda[]*a[];
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats poisson (struct Poisson p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x) {
    const vector alpha[] = {1.,1.};
    p.alpha = alpha;
  }
  if (p.lambda) {
    scalar lambda = p.lambda;
    restriction ({lambda});
  }
  else {
    const scalar lambda[] = 0.;
    p.lambda = lambda;
  }

  /**
  If the [metric](/Basilisk C#metric) is not purely Cartesian, we
  allocate and define a temporary field for the face coefficients. */

  if (constant(fm.x) != 1. || constant(fm.y) != 1.) {
    (const) face vector alpha = p.alpha;
    face vector alpha_m = new face vector;
    foreach_face()
      alpha_m.x[] = alpha.x[]*fm.x[];
    p.alpha = alpha_m;
  }
  
  /**
  We need $\alpha$ on all levels of the grid. */

  face vector alpha = p.alpha;
  restriction ((scalar *){alpha});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual, relax, &p);

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  /**
  We free the temporary face coefficients if necessary. */

  if (constant(fm.x) != 1. || constant(fm.y) != 1.)
    delete ((scalar *){alpha});

  return s;
}

/**
## Projection of a velocity field

The function below "projects" the velocity field *u* onto the space of
divergence-free velocity fields i.e.
$$
\mathbf{u}^{n+1} \leftarrow \mathbf{u} - \Delta t\alpha\nabla p
$$
so that
$$
\nabla\cdot\mathbf{u}^{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_*}{\Delta t}
$$ */

trace
mgstats project (face vector u, scalar p, (const) face vector alpha, double dt)
{

  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_*$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. *fm* is the face vector defining
  the [metric](metric.h). */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1,0] - u.x[];
    div[] /= dt*Delta;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha, tolerance = TOLERANCE/sq(dt));

  /**
  And compute $\mathbf{u}_{n+1}$ using $\mathbf{u}_*$ and $p$. If
  $\alpha$ is not defined we set it to one. */

  if (alpha.x)
    foreach_face()
      u.x[] -= dt*fm.x[]*alpha.x[]*(p[] - p[-1,0])/Delta;
  else
    foreach_face()
      u.x[] -= dt*fm.x[]*(p[] - p[-1,0])/Delta;
  boundary ((scalar *){u});

  return mgp;
}
