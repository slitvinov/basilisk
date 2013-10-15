/**
# Multigrid Poisson solvers

We want to solve Poisson equations of the general form
$$
\nabla\cdot (c\nabla a) = b
$$
This can be done efficiently using a multigrid solver. 

An important aspect of Poisson equations is that they are linear. This
property can be used to build better estimates of a solution by
successive *corrections* to an initial guess. If we define an
approximate solution $\tilde{a}$ as
$$
\tilde{a} + da = a
$$
where $a$ is the exact (unknown) solution, using the linearity of the
Laplacian we find that $da$ verifies
$$
\nabla\cdot(c\nabla da) = b - \nabla\cdot(c\nabla\tilde{a})
$$
where the right-hand-side is often called the *residual* of the
approximate solution $\tilde{a}$. 

## Homogeneous boundary conditions

To close the problem, we also need to define the boundary conditions
applicable to $da$. They will be the *homogeneous* equivalent of the
boundary conditions applied to $\tilde{a}$. The function below applies
the homogeneous boundary conditions of `a` to `da` for level `l` of
the multigrid hierarchy. */

static void boundary_homogeneous (scalar a, scalar da, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true)
      da[ghost] = a.boundary_homogeneous[b] (point, da);

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
  halo_prolongation (l, {da});
#endif
}

/**
## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
`a`, a residual `res`, a correction field `da` and a relaxation
function `relax`, we will provide an improved guess at the end of the
cycle. */

void mg_cycle (scalar a, scalar res, scalar da,
	       void (* relax) (scalar da, scalar res, int depth),
	       int nrelax, int minlevel)
{

/**
We first define the residual on all levels. */

  restriction ({res});

/**
We then proceed from the coarsest grid (`minlevel`) down to the finest grid. */

  for (int l = minlevel; l <= depth(); l++) {

/**
On the coarsest grid, we take zero as initial guess. */
    if (l == minlevel) {
      foreach_level_or_leaf (l)
	da[] = 0.;
    }

/**
On all other grids, we take as initial guess the approximate solution
on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level_or_leaf (l)
	da[] = (9.*coarse(da,0,0) + 
		3.*(coarse(da,child.x,0) + coarse(da,0,child.y)) + 
		coarse(da,child.x,child.y))/16.;

/**
We then apply homogeneous boundary conditions and do several
iterations of the relaxation function to refine the initial guess. */

    boundary_homogeneous (a, da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l);
      boundary_homogeneous (a, da, l);
    }
  }

/**
And finally we apply the resulting correction to `a`. */

  foreach()
    a[] += da[];
  boundary ({a});
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by `NITERMAX` and the
tolerance by `TOLERANCE` with the default values below.
*/

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
(and returns its maximum) as well as the relaxation function. */

mgstats mg_solve (scalar a, scalar b, 
		  double (* residual) (scalar a, scalar b, scalar res),
		  void (* relax) (scalar da, scalar res, int depth))
{
  scalar res[], da[];
  mgstats s = {0, 0., 0.};
  foreach () // fixme: need reduction(+:s.sum)
    s.sum += b[];

/**
Here we compute the initial residual field and its maximum. */

  s.maxres = residual (a, b, res);

/**
We then iterates until convergence or until `NITERMAX` is reached. Note
also that we force the solver to apply at least one cycle, even if the
initial residual is lower than `TOLERANCE`. */

  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.maxres > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, 4, 0);
    s.maxres = residual (a, b, res);
  }

/**
If we have reached the maximum number of iterations, we warn the user. */

  if (s.i == NITERMAX)
    fprintf (stderr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  sum: %g\n", 
	     NITERMAX, s.sum);
  return s;
}

/**
## Poisson equation with constant coefficients

We now apply the generic multigrid solver to the simplest case
$$
\nabla^2 a = b
$$
Using a standard 5-points discrete Laplacian operator, we easily get the
corresponding relaxation function. */

static void relax (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - sq(Delta)*b[])/4.;
}

/**
The equivalent residual function is also trivial to obtain in the case
of a Cartesian grid, however the case of the quadtree mesh requires
more careful consideration...
 */

static double residual (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = (a[] - a[-1,0])/Delta;
  boundary_normal ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta;
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + (4.*a[] - a[1,0] - a[-1,0] - a[0,1] - a[0,-1])/sq(Delta);
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

/**
Once these two functions are defined, writing the final function is trivial. */

mgstats poisson (scalar a, scalar b)
{
  return mg_solve (a, b, residual, relax);
}

/**
## Poisson equation with variable coefficients
*/

vector alpha;

static void relax_variable (scalar a, scalar b, int l)
{
  foreach_level_or_leaf (l)
    a[] = (alpha.x[1,0]*a[1,0] + alpha.x[]*a[-1,0] + 
	   alpha.y[0,1]*a[0,1] + alpha.y[]*a[0,-1] 
	   - sq(Delta)*b[])
    /(alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[]);
}

static double residual_variable (scalar a, scalar b, scalar res)
{
  double maxres = 0.;
#if QUADTREE
  /* conservative coarse/fine discretisation (2nd order) */
  vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1,0])/Delta;
  boundary_normal ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.x[] - g.x[1,0] + g.y[] - g.y[0,1])/Delta;
#else
  /* "naive" discretisation (only 1st order on quadtrees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + 
      ((alpha.x[1,0] + alpha.x[] + alpha.y[0,1] + alpha.y[])*a[]
       - alpha.x[1,0]*a[1,0] - alpha.x[]*a[-1,0] 
       - alpha.y[0,1]*a[0,1] - alpha.y[]*a[0,-1])
      /sq(Delta);
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  return maxres;
}

mgstats poisson_variable (scalar a, scalar b, vector c)
{
  alpha = c;
  restriction_staggered ({alpha});
  return mg_solve (a, b, residual_variable, relax_variable);
}
