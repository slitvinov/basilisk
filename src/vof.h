/* fixme: use unit cells centered on the origin for symmetry of
 mx + my = alpha */

static double line_alpha (double c, double nx, double ny)
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

static double line_area (double nx, double ny, double alpha)
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

static double rectangle_fraction (double c,
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

void reconstruction (const scalar c, vector n, scalar alpha)
{
  trash ({n, alpha});
  foreach()
    if (c[] > 0. && c[] < 1.) {
      // Youngs normal
      double nn = 0.;
      foreach_dimension() {
	n.x[] = (c[-1,1] + 2.*c[-1,0] + c[-1,-1] -
		 c[+1,1] - 2.*c[+1,0] - c[+1,-1]);
	nn += fabs(n.x[]);
      }
      // normalize
      if (nn > 0.)
	foreach_dimension()
	  n.x[] /= nn;
      else // this is a small fragment
	n.x[] = 1.;
      // intercept
      alpha[] = line_alpha (c[], n.x[], n.y[]);
    }
  boundary ({n, alpha});
}

static int facets (Point point,
		   scalar c, vector n, scalar alpha,
		   coord p[2])
{
  int i = 0;
  if (c[] > 0. && c[] < 1.) {
    coord o = {x, y};
    for (int s = 0; s <= 1; s++)
      foreach_dimension()
	if (fabs (n.y[]) > 1e-4) {
	  double a = (alpha[] - s*n.x[])/n.y[];
	  if (a >= 0. && a <= 1.) {
	    p[i].x = o.x - (1 - 2*s)*Delta/2.;
	    p[i++].y = o.y + Delta*(a - 0.5);
	  }
	}
    assert (i <= 2);
  }
  return i;
}

void output_facets (scalar c, FILE * fp)
{
  scalar alpha[];
  vector n[];
  reconstruction (c, n, alpha);
  coord p[2];
  foreach()
    if (facets (point, c, n, alpha, p) == 2)
      fprintf (fp, "%g %g\n%g %g\n\n", p[0].x, p[0].y, p[1].x, p[1].y);
  fflush (fp);
}

extern scalar * interfaces;
extern staggered vector u;
extern double dt;

foreach_dimension()
static void sweep_x (scalar c, scalar cc)
{
  vector n[];
  scalar alpha[], flux[];
  
  reconstruction (c, n, alpha);
  trash ({flux});
  foreach_face(x) {
    double un = u.x[]*dt/Delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double cf = (c[i,0] <= 0. || c[i,0] >= 1.) ? c[i,0] :
      rectangle_fraction (c[i,0], 
			  - s*n.x[i,0], n.y[i,0], 
			  alpha[i,0] - (s + 1.)/2.*n.x[i,0],
			  0, 0, s*un, 1);
    flux[] = u.x[]*cf;
  }
  //    boundary_normal ({flux});
  foreach()
    c[] += dt*(flux[] - flux[1,0] + cc[]*(u.x[1,0] - u.x[]))/Delta;
  boundary ({c});
}

event vof_advection (i = 1; i++)
{
  for (scalar c in interfaces) {
    scalar cc[];
    trash ({cc});
    foreach()
      cc[] = (c[] > 0.5); // eq. 19 from Weymouth & Yue
    // fixme: alternate directions
    sweep_x (c, cc);
    sweep_y (c, cc);
  }
}
