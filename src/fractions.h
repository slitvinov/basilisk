/**
# Volume fractions

These functions are used to maintain or define volume and surface
fractions either from an initial geometric definition or from an
existing volume fraction field. 

We will use basic geometric functions for square cut cells and the
"Mixed-Youngs-Centered" normal approximation of Ruben Scardovelli. */

#include "geometry.h"
#include "myc2d.h"

/**
## Coarsening and refinement of a volume fraction field 

On quadtrees, we need to define how to coarsen (i.e. "restrict") or
refine (i.e. "prolongate") interface definitions (see [geometry.h]()
for a basic explanation of how interfaces are defined).

For the normal to the interface, we don't use any interpolation from
coarse to fine i.e. we use straight "injection". */

#if QUADTREE
static double injection (Point point, scalar s)
{
  return coarse(s,0,0);
}

/**
Once we have the normal in the fine cell, we can compute the volume
fraction. */

static double fraction_prolongation (Point point, scalar c)
{

/**
If the parent cell is empty or full, we just use the same value for
the fine cell. */

  double cc = coarse(c,0,0);
  if (cc <= 0. || cc >= 1.)
    return cc;

/**
Otherwise, we reconstruct the interface in the parent cell. */

  coord n = mycs (parent, c);
  double alpha = line_alpha (cc, n.x, n.y);

/**
And compute the volume fraction in the quadrant of the coarse cell
matching the fine cell. We use symmetries to simplify the
combinations. */

  return rectangle_fraction (child.x*n.x, child.y*n.y, alpha, 
			     0., 0., 0.5, 0.5);
}

/**
The refinement function is almost exactly the same but applied to all
the children of the coarse cell. Unfortunately, we cannot just use
several applications of the prolongation function above, because of
performance. */

static void fraction_refine (Point point, scalar c)
{
  double cc = c[];
  if (cc <= 0. || cc >= 1.)
    for (int k = 0; k < 2; k++)
      for (int l = 0; l < 2; l++)
	fine(c,k,l) = cc;
  else {
    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n.x, n.y);
    for (int k = 0; k < 2; k++)
      for (int l = 0; l < 2; l++)
	fine(c,k,l) = rectangle_fraction ((2*k - 1)*n.x, (2*l - 1)*n.y, alpha,
					  0., 0., 0.5, 0.5);
  }
}

/**
Finally, we also need to prolongate the reconstructed value of
$\alpha$. This is done with the simple formula below. */

static double alpha_prolongation (Point point, scalar alpha)
{
  vector n = alpha.v;
  return 2.*coarse(alpha,0,0) - (child.x*n.x[] + child.y*n.y[])/2.;
}
#endif // QUADTREE

/**
## Computing volume fractions from a "levelset" function 

Initialising a volume fraction field representing an interface is not
trivial since it involves the numerical evaluation of surface
integrals.

Here we define a function which allows the approximation of these
surface integrals in the case of an interface defined by a "levelset"
function $\Phi$ sampled on the *vertices* of the grid.

By convention the "inside" of the interface corresponds to $\Phi > 0$.

The function takes the vertex scalar field $\Phi$ as input and fills
`c` with the volume fraction and `s` with the surface fractions
i.e. the fractions of the faces of the cell which are inside the
interface. 

![Volume and surface fractions](/src/figures/fractions.png) */

void fractions (const scalar Phi, scalar c, staggered vector s)
{

/**
### Surface fraction computation

We start by computing the surface fractions. */

  foreach_face() {

/**
If the values of $\Phi$ on the vertices of the face have opposite
signs, we know that the face is cut by the interface. */

    if (Phi[]*Phi[0,1] < 0.) {

/**
In that case we can find an approximation of the interface position by
simple linear interpolation. We also check the sign of one of the
vertices to orient the interface properly. */

      s.x[] = Phi[]/(Phi[] - Phi[0,1]);
      if (Phi[] < 0.)
	s.x[] = 1. - s.x[];
    }
    
/**
If the values of $\Phi$ on the vertices of the face have the same sign
(or are zero), then the face is either entirely outside or entirely
inside the interface. We check the sign of both vertices to treat
limit cases properly (when the interface intersects the face exactly
on one of the vertices). */

    else
      s.x[] = (Phi[] > 0. || Phi[0,1] > 0.);
  }

/**
We make sure that surface fractions are defined conservatively (on
quadtree meshes). */

  boundary_normal ({s});

/**
### Volume fraction computation */

  foreach() {

/**

We first compute the normal to the interface. This can be done easily
using the surface fractions. The idea is to compute the circulation of
the normal along the boundary $\partial\Omega$ of the fraction of the
cell $\Omega$ inside the interface. Since this is a closed curve, we
have
$$
\oint_{\partial\Omega}\mathbf{n}\;dl = 0
$$ 
We can further decompose the integral into its parts along the faces
of the square and the part along the interface. For the case pictured
above, we get for one component (and similarly for the other)
$$
- s_x[] + \oint_{\Phi=0}n_x\;dl = 0
$$
If we now define the *average normal* to the interface as
$$
\overline{\mathbf{n}} = \oint_{\Phi=0}\mathbf{n}\;dl
$$
We have in the general case
$$
\overline{\mathbf{n}}_x = s_x[] - s_x[1,0]
$$
and
$$
|\overline{\mathbf{n}}| = \oint_{\Phi=0}\;dl
$$ 
Note also that this average normal is exact in the case of a linear
interface. */

    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1,0];
      nn += fabs(n.x);
    }

/**
If the norm is zero, the cell is full or empty and the volume fraction
is identical to one of the surface fractions. */

    if (nn == 0.)
      c[] = s.x[];

/**
Otherwise we are in a cell containing the interface. We first
normalise the normal. */

    else {
      foreach_dimension()
	n.x /= nn;

/**
To find the intercept $\alpha$, we look for a face which is cut by the
interface, find the coordinate $a$ of the intersection and use it to
derive $\alpha$. */

      double alpha = undefined;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (s.x[i,0] > 0. && s.x[i,0] < 1.) {
	    double a = sign(Phi[i,0])*(s.x[i,0] - 0.5);
	    alpha = n.x*(i - 0.5) + n.y*a;
	  }

/**
Once we have $\mathbf{n}$ and $\alpha$, the (linear) interface is
fully defined and we can compute the volume fraction using our
pre-defined function. */

      c[] = line_area (n.x, n.y, alpha);
    }
  }

/**
On a quadtree grid, we set the prolongation and refinement functions
we defined above. */

#if QUADTREE
  c.prolongation = fraction_prolongation;
  c.refine = fraction_refine;
#endif

/**
Finally we apply the boundary conditions. */

  boundary ({c});
}

/**
## Interface reconstruction from volume fractions

The reconstruction of the interface geometry from the volume fraction
field requires computing an approximation to the interface normal.

### Youngs normal approximation 

This a simple, but relatively inaccurate way of approximating the
normal. It is simply a weighted average of centered volume fraction
gradients. We include it as an example but it is not used. */

coord youngs_normal (Point point, scalar c)
{
  coord n;
  double nn = 0.;
  foreach_dimension() {
    n.x = (c[-1,1] + 2.*c[-1,0] + c[-1,-1] -
	   c[+1,1] - 2.*c[+1,0] - c[+1,-1]);
    nn += fabs(n.x);
  }
  // normalize
  if (nn > 0.)
    foreach_dimension()
      n.x /= nn;
  else // this is a small fragment
    n.x = 1.;
  return n;
}

/**
### Interface reconstruction 

The reconstruction function takes a volume fraction field `c` and
returns the corresponding normal vector field `n` and intercept field
$\alpha$. */

void reconstruction (const scalar c, vector n, scalar alpha)
{
  foreach() {

/**
If the cell is empty or full, we set $\mathbf{n}$ and $\alpha$ only to
avoid using unitialised values in `alpha_prolongation()`. */

    if (c[] <= 0. || c[] >= 1.)
      alpha[] = n.x[] = n.y[] = 0.;
    else {

/**
Otherwise, we compute the interface normal using the
Mixed-Youngs-Centered scheme, copy the result into the normal field
and compute the intercept $\alpha$ using our predefined function. */

      coord m = mycs (point, c);
      // coord m = youngs_normal (point, c);
      foreach_dimension()
	n.x[] = m.x;
      alpha[] = line_alpha (c[], m.x, m.y);
    }
  }

#if QUADTREE

/**
On a quadtree grid, we set the prolongation functions defined above.
We do not restrict the normal or the intercept (they are recomputed
from the volume fraction when needed on the coarse mesh) */

  n.x.coarsen = n.y.coarsen = alpha.coarsen = coarsen_none;
  n.x.prolongation = n.y.prolongation = injection;

/**
For $\alpha$ we store the normal field in the `v` attribute (which is
not really meant to be used this way) to be able to pass it to the
prolongation function. */

  alpha.v = n;
  alpha.prolongation = alpha_prolongation;
#endif

/**
Finally we apply the boundary conditions to define $\mathbf{n}$ and
$\alpha$ everywhere (using the prolongation functions when necessary
on quadtree grids). */

  boundary ({n, alpha});
}

/**
## Interface outputs 

This first function takes a volume (resp. face) fraction field `c`
(resp. `s`) and "draws" the corresponding interfaces facets in file
`fp`. The segment endpoints are defined by pairs of coordinates. Each
pair of endpoints is separated from the next pair by a newline, so
that the resulting file is directly visualisable with gnuplot.

Note that this function will produce a continuous interface representation 
in most cases. */

void output_fractions (const scalar c, const staggered vector s, FILE * fp)
{
  foreach() {
    // compute normal from fractions
    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1,0];
      nn += fabs(n.x);
    }
    if (nn > 0.) { // interfacial cell
      foreach_dimension()
	n.x /= nn;
      double alpha = line_alpha (c[], n.x, n.y);
      coord segment[2];
      if (facets (c[], n, alpha, segment) == 2)
	fprintf (fp, "%g %g\n%g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 x + segment[1].x*Delta, y + segment[1].y*Delta);
    }
  }
  fflush (fp);
}

/**
This second function only takes a volume fraction field and uses
interface reconstruction to derive the full geometry necessary to draw
the facets.

In most cases it will only produce a piecewise continuous (i.e. VOF)
interface representation. */

void output_facets (scalar c, FILE * fp)
{
  scalar alpha[];
  vector n[];
  reconstruction (c, n, alpha);
  coord segment[2];
  foreach() {
    coord m = {n.x[],n.y[]};
    if (facets (c[], m, alpha[], segment) == 2)
      fprintf (fp, "%g %g\n%g %g\n\n", 
	       x + segment[0].x*Delta, y + segment[0].y*Delta, 
	       x + segment[1].x*Delta, y + segment[1].y*Delta);
  }
  fflush (fp);
}
