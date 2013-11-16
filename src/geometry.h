/* fixme: use unit cells centered on the origin for symmetry of
 mx + my = alpha */

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

  return alpha;
}

double line_area (double nx, double ny, double alpha)
{
  double alpha1, a, v, area;

  alpha1 = alpha;
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

typedef struct {
  double x, y;
} coord;

double rectangle_fraction (double c,
			   double nx, double ny, double alpha,
			   double x1, double y1,
			   double x2, double y2)
{
  coord n = {nx, ny}, r[2] = {{x1,y1},{x2,y2}};
  foreach_dimension() {
    alpha -= n.x*r[0].x;
    n.x *= r[1].x - r[0].x;
  }
  return line_area (n.x, n.y, alpha);
}

int facets (double c, coord n, double alpha,
	    coord p[2])
{
  int i = 0;
  if (c > 0. && c < 1.) {
    for (int s = 0; s <= 1; s++)
      foreach_dimension()
	if (fabs (n.y) > 1e-4) {
	  double a = (alpha - s*n.x)/n.y;
	  if (a >= 0. && a <= 1.) {
	    p[i].x   = s - 0.5;
	    p[i++].y = a - 0.5;
	  }
	}
    assert (i <= 2);
  }
  return i;
}
