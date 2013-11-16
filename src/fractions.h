#include "geometry.h"

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
      n.x /= nn; n.y /= nn;
      double alpha = undefined;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (s.x[i,0] > 0. && s.x[i,0] < 1.) {
	    double a = phi[i,0] > 0. ? s.x[i,0] : 1. - s.x[i,0];
	    alpha = n.x*i + n.y*a;
	  }
      c[] = line_area (n.x, n.y, alpha);
    }
  }
  boundary ({c});
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
      n.x /= nn; n.y /= nn;
      double alpha = line_alpha (c[], n.x, n.y);
      coord p[2];
      if (facets (c[], n, alpha, p) == 2)
	fprintf (fp, "%g %g\n%g %g\n\n", 
		 x + p[0].x*Delta, y + p[0].y*Delta, 
		 x + p[1].x*Delta, y + p[1].y*Delta);
    }
  }
  fflush (fp);
}

void reconstruction (const scalar c, vector n, scalar alpha)
{
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

void output_facets (scalar c, FILE * fp)
{
  scalar alpha[];
  vector n[];
  reconstruction (c, n, alpha);
  coord p[2];
  foreach() {
    coord m = {n.x[],n.y[]};
    if (facets (c[], m, alpha[], p) == 2)
      fprintf (fp, "%g %g\n%g %g\n\n", 
	       x + p[0].x*Delta, y + p[0].y*Delta, 
	       x + p[1].x*Delta, y + p[1].y*Delta);
  }
  fflush (fp);
}
