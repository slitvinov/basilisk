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

## Homogeneous boundary conditions

To close the problem, we also need to define the boundary conditions
applicable to $da$. They will be the *homogeneous* equivalent of the
boundary conditions applied to $\tilde{a}$. The function below applies
the homogeneous boundary conditions of `a` to `da` for level `l` of
the multigrid hierarchy. */

static void boundary_homogeneous (scalar * a, scalar * da, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) {
      scalar s, ds;
      for (s, ds in a, da)
	ds[ghost] = s.boundary_homogeneous[b] (point, ds);
    }

/**
On quadtree meshes, we also need to make sure that stencils are
consistent. In the general case, this involves first restricting the
field values from the fine mesh to the coarse mesh and then
prolongating the field on 'halo' cells using the coarse grid
stencils. In the case of the multigrid cycle below, the solution is
already known on coarse grids so that the restriction operation is not
necessary, which leaves only the prolongation operation for level
`l`. */

#if QUADTREE
  halo_prolongation (l, da);
#endif
}

/**
## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
`a`, a residual `res`, a correction field `da` and a relaxation
function `relax`, we will provide an improved guess at the end of the
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
We then proceed from the coarsest grid (`minlevel`) down to the finest grid. */

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
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  s[] = (9.*coarse(s,0,0) + 
		 3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;

/**
We then apply homogeneous boundary conditions and do several
iterations of the relaxation function to refine the initial guess. */

    boundary_homogeneous (a, da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_homogeneous (a, da, l);
    }
  }

/**
And finally we apply the resulting correction to `a`. */

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

The maximum number of iterations is controlled by `NITERMAX` and the
tolerance by `TOLERANCE` with the default values below. */

int NITERMAX = 100;
double TOLERANCE = 1e-3;

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double maxres, sum; // maximum residual, sum of r.h.s.
} mgstats;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer `data` can be used to pass arguments to these
functions. */

mgstats mg_solve (scalar * a, scalar * b,
		  double (* residual) (scalar * a, scalar * b, scalar ** res,
				       void * data),
		  void (* relax) (scalar * da, scalar * res, int depth, 
				  void * data),
		  void * data)
{

/**
The list of residual fields will be allocated by the residual()
function. We allocate a new correction field for each of the scalars
in `a`. */

  scalar * res = NULL, * da = NULL;
  for (scalar s in a) {
    scalar ds = new scalar;
    da = list_append (da, ds);
  }

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

  s.maxres = residual (a, b, &res, data);

/**
We then iterates until convergence or until `NITERMAX` is reached. Note
also that we force the solver to apply at least one cycle, even if the
initial residual is lower than `TOLERANCE`. */

  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.maxres > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, data, 4, 0);
    s.maxres = residual (a, b, &res, data);
  }

/**
If we have reached the maximum number of iterations, we warn the user. */

  if (s.i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  sum: %g\n", 
	     NITERMAX, s.sum);

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

`alpha` and `lambda` are declared as `(const)` to indicate that the
function works also when `alpha` and `lambda` are constant vector
(resp. scalar) fields. */

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
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
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the quadtree mesh
requires more careful consideration... */

static double residual (scalar * al, scalar * bl, scalar ** resl, void * data)
{
  scalar a = al[0], b = bl[0];

/**
If the residual field is not already allocated we allocate it as well
as the associated list. */

  scalar res;
  if (*resl)
    res = (*resl)[0];
  else {
    res = new scalar;
    *resl = list_append (NULL, res);
  }

  struct Poisson * p = data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1,0])/Delta;
  boundary_normal ({g});
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

  if (p.alpha.x) {
    face vector alpha = p.alpha;
    restriction ((scalar *){alpha});
  }
  else {
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
  
  scalar a = p.a, b = p.b;
  return mg_solve ({a}, {b}, residual, relax, &p);
}
