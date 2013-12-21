/**
# A generic solver for systems of conservation laws

Using the ideas of [Kurganov and Tadmor,
2000](references.bib#kurganov2000) it is possible to write a generic
solver for systems of conservation laws of the form
$$
\partial_t\left(\begin{array}{c}
    s_i\\
    \mathbf{v}_j\\
 \end{array}\right) + \nabla\cdot\left(\begin{array}{c}
    \mathbf{F}_i\\
    \mathbf{T}_j\\
 \end{array}\right) = 0
$$
where $s_i$ is a list of scalar fields, $\mathbf{v}_j$ a list of
vector fields and $\mathbf{F}_i$, $\mathbf{T}_j$ are the corresponding
vector (resp. tensor) fluxes. 

Note that the [Saint-Venant solver](saint-venant.h) is a particular
case of this generic algorithm.

The user must provide the lists of conserved scalar and vector fields
*/

extern scalar * scalars;
extern vector * vectors;

/**
as well as a function which, given the state of each quantity,
returns the fluxes and the minimum/maximum eigenvalues
(i.e. `eigenvalue[0]`/`eigenvalue[1]`) of the corresponding Riemann
problem. */

void flux (const double * state, double * flux, double * eigenvalue);

/**
Optional source terms can be added by overloading the *sources* function
pointer. */

static void no_sources (scalar * current, scalar * updates)
{
  foreach()
    for (scalar s in updates)
      s[] = 0.;
}

void (* sources) (scalar * current, scalar * updates) = no_sources;

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in `predictor-corrector.h` needs
to know which fields are updated i.e. all the scalars + the components
of all the vector fields. */

scalar * evolving;

event defaults (i = 0)
{
  evolving = list_concat (scalars, (scalar *) vectors);
}

/**
We force boundary conditions on all fields after initialisation.
*/

event init (i = 0)
{
  boundary (all);
}

/**
### Computing fluxes

The core of the central-upwind scheme (see e.g. section 3.1 of
[Kurganov & Levy, 2002](references.bib#kurganov2002)) is the
approximate solution of the Riemann problem given by the left and
right states to get the fluxes `f`. */

static double riemann (const double * right, const double * left,
		       double Delta, double * f, int len, 
		       double dtmax)
{
  double fm[len], fp[len], em[2], ep[2];
  flux (right, fm, em);
  flux (left,  fp, ep);
  double ap = max(ep[1], em[1]); ap = max(ap, 0.);
  double am = min(ep[0], em[0]); am = min(am, 0.);
  double a = max(ap, -am);
  if (a > 0.) {
    for (int i = 0; i < len; i++)
      f[i] = (ap*fm[i] - am*fp[i] + ap*am*(left[i] - right[i]))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < dtmax)
      dtmax = dt;
  }
  else
    for (int i = 0; i < len; i++)
      f[i] = 0.;
  return dtmax;
}

double update (scalar * conserved, scalar * updates, double dtmax)
{

/**
We declare an empty list of fluxes which will be filled with fluxes
for each scalar and vector quantity. */

  vector * lflux = NULL;

/**
The gradients of each quantity are stored in a list of dynamically-allocated
fields. First-order reconstruction is used for the gradient fields. */

  vector * slopes = NULL;
  for (scalar s in conserved) {
    vector slope = new vector;
    foreach_dimension()
      slope.x.gradient = zero;
    slopes = vectors_append (slopes, slope);
  }
  gradients (conserved, slopes);

/**
The flux fields are dynamically allocated and the left and right
states of each conserved quantity are locally-allocated. */

  int len = list_len (conserved);
  double r[len], l[len]; // right/left Riemann states
  double f[len];         // fluxes for each conserved quantity
  for (scalar s in conserved) {
    vector f1 = new vector;
    lflux = vectors_append (lflux, f1);
  }

/**
The predictor-corrector scheme treats all fields as scalars (stored in
the `conserved` list). We need to recover vector and tensor quantities
from these lists. To do so, knowing the number of scalar fields, we
split the scalar list into a list of scalars and a list of vectors. */

  int scalars_len = list_len (scalars);

  scalar * scalars = list_copy (conserved);
  if (scalars) scalars[scalars_len] = -1;
  vector * vectors = vectors_from_scalars (&conserved[scalars_len]);

/**
We then do the same for the gradients i.e. split the list of vectors
into a list of vectors and a list of tensors. */

  vector * scalar_slopes = vectors_copy (slopes);
  if (scalar_slopes) scalar_slopes[scalars_len] = (vector){-1,-1};
  tensor * vector_slopes = tensors_from_vectors (&slopes[scalars_len]);

/**
And again for the fluxes. */

  vector * scalar_fluxes = vectors_copy (lflux);
  if (scalar_fluxes) scalar_fluxes[scalars_len] = (vector){-1,-1};
  tensor * vector_fluxes = tensors_from_vectors (&lflux[scalars_len]);

/**
We are ready to compute the fluxes through each face of the domain. */

  foreach_face (reduction (min:dtmax)) {

/**
#### Left/right state reconstruction 

We use the central values of each scalar/vector quantity and the
pre-computed gradients to compute the left and right states. */

    double dx = Delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in scalars,scalar_slopes) {
      l[i] = s[] - dx*g.x[];
      r[i++] = s[-1,0] + dx*g.x[-1,0];
    }
    vector v;
    tensor t;
    for (v,t in vectors,vector_slopes) {
      l[i] = v.x[] - dx*t.x.x[];
      r[i++] = v.x[-1,0] + dx*t.x.x[-1,0];
      l[i] = v.y[] - dx*t.y.x[];
      r[i++] = v.y[-1,0] + dx*t.y.x[-1,0];
    }

/**
#### Riemann problem

We then call the generic Riemann solver and store the resulting fluxes
in the pre-allocated fields. */

    dtmax = riemann (r, l, Delta, f, len, dtmax);
    i = 0;
    for (vector fs in scalar_fluxes)
      fs.x[] = f[i++];
    for (tensor fv in vector_fluxes) {
      fv.x.x[] = f[i++];
      fv.y.x[] = f[i++];
    }
  }

  boundary_normal (lflux);

/**
#### Update

The update for each scalar quantity is the divergence of the fluxes. */

  sources (conserved, updates);
  foreach() {
    scalar ds;
    vector f;
    for (ds,f in updates,lflux)
      ds[] += (f.x[] - f.x[1,0] + f.y[] - f.y[0,1])/Delta;
  }

/**
#### Cleanup

We finally deallocate the memory used to store lists and gradient
fields. */

  free (scalars);
  free (vectors);
  free (scalar_slopes);
  free (vector_slopes);
  free (scalar_fluxes);
  free (vector_fluxes);
  delete ((scalar *) slopes);
  free (slopes);
  delete ((scalar *) lflux);
  free (lflux);

  return dtmax;
}
