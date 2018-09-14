/**
# Hydrostatic balance with refined solids
*/

#include "solid.h"
#include "navier-stokes/centered.h"

void circle (scalar cs, face vector fs, double xs, double ys)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = - (sq(x) + sq(y) - sq(0.3));
  fractions (phi, cs, fs);
  restriction ({cs, fs});
}

void porous (scalar cs, face vector fs, double xs, double ys)
{
#if 1
  int ns = 200; // 160, 80
  double xc[ns], yc[ns], R[ns];
  srand (0);
  for (int i = 0; i < ns; i++)
    xc[i] = xs + 0.5*noise(),
      yc[i] = ys + 0.5*noise(), R[i] = 0.02 + 0.04*fabs(noise());
#else
  int ns = 1;
  double xc[ns], yc[ns], R[ns];
  xc[0] = xs; yc[0] = ys; R[0] = 0.3;
#endif
  
  vertex scalar phi[];
  int na = 0;
  do {
    foreach_vertex() {
      phi[] = HUGE;
      for (double xp = -L0; xp <= L0; xp += L0)
	for (double yp = -L0; yp <= L0; yp += L0)
	  for (int i = 0; i < ns; i++)
	    phi[] = intersection (phi[], (sq(x + xp - xc[i]) +
					  sq(y + yp - yc[i]) - sq(R[i])));
    }
    fractions (phi, cs, fs);
    na++;
  } while (0); // adapt_wavelet ({cs}, (double[]){1e-2}, 8).nf && na < 10);

  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -L0; xp <= L0; xp += L0)
      for (double yp = -L0; yp <= L0; yp += L0)
	for (int i = 0; i < ns; i++)
	  phi[] = intersection (phi[], (sq(x + xp - xc[i]) +
					sq(y + yp - yc[i]) - sq(R[i])));
  }
  fractions (phi, cs, fs);

#if 1
  foreach_face()
#if 0
    if (sq(fs.x[]) > 20.*cs[] || sq(fs.x[]) > 20.*cs[-1])
      fs.x[] = 0.;
#else
    if (!cs[] || !cs[-1])
      fs.x[] = 0.;
#endif
  boundary ((scalar *){cs, fs});
#endif
  restriction ({cs, fs});

  for (int i = 0; i <= depth(); i++)
    foreach_level (i)
      if (!cs[])
	assert (!fs.x[] && !fs.x[1] && !fs.y[] && !fs.y[0,1]);
}

static inline void restriction_solid_linear (Point point, scalar s)
{
  // 0 children
  if (!cs[])
    return;
  // 4 children
  if (fine(cs,0,0) && fine(cs,1,0) && fine(cs,0,1) && fine(cs,1,1)) {
    s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
    return;
  }
  // 3 children
  if (fine(cs,0,0) && fine(cs,1,1)) {
    s[] = (fine(s,0,0) + fine(s,1,1))/2.;
    return;
  }
  if (fine(cs,0,1) && fine(cs,1,0)) {
    s[] = (fine(s,0,1) + fine(s,1,0))/2.;
    return;
  }
  // 2 children
  foreach_dimension()
    for (int i = 0; i <= 1; i++)
      if (fine(cs,i,0) && fine(cs,i,1)) {
	if (is_leaf(neighbor(2*i-1))) {
	  assert (cs[2*i-1]);
	  s[] = (2.*(fine(s,i,0) + fine(s,i,1)) - s[2*i-1])/3.;
	}
	else {
	  double v;
	  if (fine(cs,3*i-1,0)) {
	    if (fine(cs,3*i-1,1))
	      v = fine(s,3*i-1,0) + fine(s,3*i-1,1);
	    else
	      v = 2.*fine(s,3*i-1,0);
	  }
	  else if (fine(cs,3*i-1,1))
	    v = 2.*fine(s,3*i-1,1);
	  else {
	    s[] = (fine(s,i,0) + fine(s,i,1))/2.;
	    return;
	  }	    
	  s[] = (3.*(fine(s,i,0) + fine(s,i,1)) - v)/4.;
	}
	return;
      }
  // 1 child
  for (int i = 0; i <= 1; i++)
    for (int j = 0; j <= 1; j++)
      if (fine(cs,i,j)) {
	if (is_leaf(neighbor(2*i-1,2*j-1))) {
	  assert (cs[2*i-1,2*j-1]);
	  s[] = (4.*fine(s,i,j) - s[2*i-1,2*j-1])/3.;
	}
	else if (fine(cs,3*i-1,3*j-1))
	  s[] = (3.*fine(s,i,j) - fine(s,3*i-1,3*j-1))/2.;
	else
	  s[] = fine(s,i,j);
	return;
      }
  assert (false);
}

static inline void refine_solid_linear (Point point, scalar s)
{
  if (!cs[])
    return;
  foreach_child()
    if (cs[]) {
      if (coarse(cs,child.x) && coarse(cs,0,child.y) &&
	  coarse(cs,child.x,child.y))
	// bilinear interpolation
	s[] = (9.*coarse(s) + 
	       3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
	       coarse(s,child.x,child.y))/16.;
      else if (coarse(cs,child.x) && coarse(cs,0,child.y))
	// triangular interpolation
	s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      else {
	// pathological cases
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(cs,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (!is_leaf(aparent(0,child.y)) &&
		   aparent(0,child.y).neighbors &&
		   cs[0,child.y] && cs[-child.x,child.y])
	    s[] += (s[0,child.y] - s[-child.x,child.y])/2.;
	}
      }
    }
}

int main()
{
  origin (-0.5, -0.5);

  init_grid (1 << 6);

  event ("defaults");
  event ("metric");

  porous (cs, fs, 0, 0);
#if 1
  refine (level < 8 && cs[] > 0 && cs[] < 1);
  porous (cs, fs, 0, 0);
#endif

  const face vector G[] = {1.,1.};
  a = G;
  
  TOLERANCE = 1e-6;
  NITERMAX = 100;
  mgp.nrelax = 100;
  alpha = fm;
  dt = 1.;

  p.restriction = restriction_solid_linear;
  p.prolongation = p.refine = refine_solid_linear;

  event ("acceleration");
#if 1
  event ("projection");
#else
  foreach()
    p[] = G.x[]*x + G.y[]*y; // exact pressure
  boundary ({p});
  foreach_face()
    uf.x[] -= alpha.x[] ? dt*alpha.x[]*face_gradient_x (p, 0) : 0.;
  boundary ((scalar *){uf});
  
  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[] ? fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta : 0.;
  boundary_flux ({gf});
  
  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
  boundary ((scalar *){g});

  correction (dt);
#endif
  
  fprintf (stderr, "mgp %g %g %d %d\n", mgp.resb, mgp.resa, mgp.i, mgp.minlevel);
  fprintf (stderr, "umax %g %g\n", normf(u.x).max, normf(u.y).max);

  scalar p1[];
  vector u1[];
  foreach()
    p1[] = p[], u1.x[] = uf.x[], u1.y[] = uf.y[];
  dump();
}
