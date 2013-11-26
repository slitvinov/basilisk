/**
# Basic geometric functions

These basic geometric functions are mostly related to Volume-Of-Fluid
computations.

We consider a square cell of size unity centered on the origin, cut by
a straight line.

![Cell and interface](/src/figures/square.png)

The line can be described by the equation
$$
n_xx+n_yy=\alpha
$$
where $\mathbf{n}$ is a vector normal to the interface and $\alpha$ is
the intercept. We note $c$ the volume of the part of the square cell
which lies "inside" the interface, where "inside" is defined by
convention as the opposite direction to the normal vector $\mathbf{n}$
(i.e. the normal vector is pointing "outside").

With these definitions, the interface is uniquely defined by providing
$\mathbf{n}$ and either $\alpha$ or $c$ i.e. there is a unique
function which computes $\alpha$ given $c$ and $\mathbf{n}$. We call
this function `line_alpha()` and define it as: */

double line_alpha (double c, double nx, double ny)
{
  double alpha, n1, n2;
  
  n1 = fabs (nx); n2 = fabs (ny);
  if (n1 > n2)
    swap (double, n1, n2);
  
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (nx < 0.)
    alpha += nx;
  if (ny < 0.)
    alpha += ny;

  return alpha - (nx + ny)/2.;
}

/**
Conversely there is a unique function computing $c$ as a function of
$\mathbf{n}$ and $\alpha$. We call this function `line_area()` and
define it as: */

double line_area (double nx, double ny, double alpha)
{
  double alpha1, a, v, area;

  alpha1 = alpha + (nx + ny)/2.;
  if (nx < 0.) {
    alpha1 -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha1 -= ny;
    ny = - ny;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= nx + ny)
    return 1.;

  if (nx == 0.)
    area = alpha1/ny;
  else if (ny == 0.)
    area = alpha1/nx;
  else {
    v = alpha1*alpha1;

    a = alpha1 - nx;
    if (a > 0.)
      v -= a*a;
    
    a = alpha1 - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}

/**
VOF algorithms require the computation of volume fractions on
(rectangular) parts of the initial square cell.

We first define a new type used to store the coordinates of the
vertices of the rectangle. */

typedef struct {
  double x, y;
} coord;

/**
We then define a function which takes an interface definition
($\mathbf{n}$, $\alpha$), the coordinates of the lower-left `(x1,y1)`
and upper-right `(x2,y2)` corners of a rectangle and returns the
fraction of this rectangle which lies inside the interface. */

double rectangle_fraction (double nx, double ny, double alpha,
			   double x1, double y1,
			   double x2, double y2)
{
  coord n = {nx, ny}, r[2] = {{x1,y1},{x2,y2}};
  foreach_dimension() {
    alpha -= n.x*(r[1].x + r[0].x)/2.;
    n.x *= r[1].x - r[0].x;
  }
  return line_area (n.x, n.y, alpha);
}

/**
From the interface definition, it is also possible to compute the
coordinates of the segment representing the interface in the unit
cell.

The function below returns the 0,1 or 2 coordinates (stored in the `p`
array provided by the user) of the corresponding interface
segments. The case where only 1 coordinate is returned corresponds to
the degenerate case where the interface intersects the cell exactly on
a vertex. */

int facets (double c, coord n, double alpha,
	    coord p[2])
{
  int i = 0;
  if (c > 0. && c < 1.) {
    for (double s = -0.5; s <= 0.5; s += 1.)
      foreach_dimension()
	if (fabs (n.y) > 1e-4) {
	  double a = (alpha - s*n.x)/n.y;
	  if (a >= -0.5 && a <= 0.5) {
	    p[i].x   = s;
	    p[i++].y = a;
	  }
	}
    assert (i <= 2);
  }
  return i;
}
