#include "geometry.h"
#include "myc2d.h"

#if QUADTREE
static double injection (Point point, scalar s)
{
  return coarse(s,0,0);
}

static double fraction_prolongation (Point point, scalar c)
{
  double cc = coarse(c,0,0);
  if (cc <= 0. || cc >= 1.)
    return cc;
  coord n = mycs (parent, c);
  double alpha = line_alpha (cc, n.x, n.y);
  return rectangle_fraction (child.x*n.x, child.y*n.y, alpha, 
			     0., 0., 0.5, 0.5);
}

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

static double alpha_prolongation (Point point, scalar alpha)
{
  vector n = alpha.v;
  return 2.*coarse(alpha,0,0) - (child.x*n.x[] + child.y*n.y[])/2.;
}
#endif // QUADTREE

void fractions (const scalar phi, scalar c, staggered vector s)
{
  foreach_face() {
    if (phi[]*phi[0,1] < 0.) {
      s.x[] = phi[]/(phi[] - phi[0,1]);
      if (phi[] < 0.)
	s.x[] = 1. - s.x[];
    }
    else
      s.x[] = (phi[] > 0. || phi[0,1] > 0.);
  }
  boundary_normal ({s});
  foreach() {
    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1,0];
      nn += fabs(n.x);
    }
    if (nn == 0.) // full or empty cell
      c[] = s.x[];
    else { // interfacial cell
      foreach_dimension()
	n.x /= nn;
      double alpha = undefined;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (s.x[i,0] > 0. && s.x[i,0] < 1.) {
	    double a = sign(phi[i,0])*(s.x[i,0] - 0.5);
	    alpha = n.x*(i - 0.5) + n.y*a;
	  }
      c[] = line_area (n.x, n.y, alpha);
    }
  }
#if QUADTREE
  c.prolongation = fraction_prolongation;
  c.refine = fraction_refine;
#endif
  boundary ({c});
}

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

void reconstruction (const scalar c, vector n, scalar alpha)
{
  foreach() {
    if (c[] <= 0. || c[] >= 1.)
      // this is just so that alpha_prolongation() does not blow up
      alpha[] = n.x[] = n.y[] = 0.;
    else {
      coord m = mycs (point, c); // mixed Youngs/centered normal
      // coord m = youngs_normal (point, c);
      foreach_dimension()
	n.x[] = m.x;
      // intercept
      alpha[] = line_alpha (c[], m.x, m.y);
    }
  }
#if QUADTREE
  n.x.coarsen = n.y.coarsen = alpha.coarsen = coarsen_none;
  n.x.prolongation = n.y.prolongation = injection;
  alpha.v = n;
  alpha.prolongation = alpha_prolongation;
#endif
  boundary ({n, alpha});
}

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
