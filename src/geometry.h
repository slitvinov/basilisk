/**
# Basic geometric functions

These basic geometric functions are mostly related to Volume-Of-Fluid
computations.

We consider a square cell of size unity centered on the origin, cut by
a straight line.

![Cell and interface](/src/figures/square.svg)

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

double line_alpha (double c, coord n)
{
  double alpha, n1, n2;
  
  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);
  
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
#else // dimension == 3
  double alpha;
  coord n1;
  
  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p = 2.*m1*m2;
    double q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    double p12 = sqrt (p);
    double teta = acos(q/(p*p12))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4.;
    double q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    double p12 = sqrt(p);
    double teta = acos(q/(p*p12))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha;
#endif // dimension == 3
}

/**
Conversely there is a unique function computing $c$ as a function of
$\mathbf{n}$ and $\alpha$. We call this function `line_area()` and
define it as: */

double line_area (coord n, double alpha)
{
  double alpha1, a, v, area;

  alpha1 = alpha + (n.x + n.y)/2.;
  if (n.x < 0.) {
    alpha1 -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha1 -= n.y;
    n.y = - n.y;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y)
    return 1.;

  if (n.x == 0.)
    area = alpha1/n.y;
  else if (n.y == 0.)
    area = alpha1/n.x;
  else {
    v = alpha1*alpha1;

    a = alpha1 - n.x;
    if (a > 0.)
      v -= a*a;
    
    a = alpha1 - n.y;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*n.x*n.y);
  }

  return clamp (area, 0., 1.);
#else // dimension == 3
  double al = alpha + max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  g_assert (tmp > 0.);
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1*1, n2);
  double b3 = max(n1*1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
	   b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) + 
	   b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
#endif // dimension == 3
}

/**
VOF algorithms require the computation of volume fractions on
(rectangular) parts of the initial square cell.

We first define a function which takes an interface definition
($\mathbf{n}$, $\alpha$), the coordinates of the lower-left `(a_x,a_y)`
and upper-right `(b_x,b_y)` corners of a rectangle and returns the
fraction of this rectangle which lies inside the interface. */

double rectangle_fraction (double n_x, double n_y, double n_z,
			   double alpha,
			   double a_x, double a_y, double a_z,
			   double b_x, double b_y, double b_z)
{
  coord n;
  foreach_dimension() {
    alpha -= n_x*(b_x + a_x)/2.;
    n.x = n_x*(b_x - a_x);
  }
  return line_area (n, alpha);
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
	if (fabs (n.y) > 1e-4 && i < 2) {
	  double a = (alpha - s*n.x)/n.y;
	  if (a >= -0.5 && a <= 0.5) {
	    p[i].x   = s;
	    p[i++].y = a;
	  }
	}
  }
  return i;
}
