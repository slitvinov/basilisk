#define NH 1
#define BOX 1 // fixme

#include "poisson.h"

scalar * qzl = NULL, * phil = NULL;
mgstats mgp;

event defaults (i = 0)
{
  non_hydro = true;
  
  assert (qzl == NULL && phil == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar qz = new scalar;
    scalar phi = new scalar;
    qzl = list_append (qzl, qz);
    phil = list_append (phil, phi);
  }
  reset (qzl, 0.);
  reset (phil, 0.);

  if (!tracers)
    tracers = calloc (nl, sizeof(scalar *));

  int l = 0;
  for (scalar qz in qzl) {
    tracers[l] = list_append (tracers[l], qz);
    l++;
  }
}

event viscous_term (i++)
{
  if (nu > 0.) {
    // fixme: ugly hack
    scalar lb = lambda_b, d = dut, u = u_b;
    lambda_b = dut = u_b = zeroc;
    foreach()
      // fixme: BCs should be different from those of horizontal velocity
      vertical_viscosity (point, hl, (scalar *) qzl, dt);
    boundary (qzl);
    lambda_b = lb; dut = d; u_b = u;
  }
}

static void box_matrix (Point point, scalar * phil, scalar * rhsl,
			double * H, double * d)
{
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    scalar h = hl[m], rhs = rhsl[m], phi = phil[m];
    d[l] = rhs[];
    foreach_dimension()
      d[l] -= h[]*(h[-1]*phi[-1] + h[1]*phi[1])/sq(Delta);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    H[l*nl + l] = - 2.*dimension*sq(h[])/sq(Delta) - 4.;
    if (l > 0) {
      scalar phip = phil[m+1];
      foreach_dimension()
	d[l] -= h[]*(h[-1]*phip[-1] + h[1]*phip[1])/sq(Delta);
      H[l*(nl + 1) - 1] = - 2.*dimension*sq(h[])/sq(Delta) + 4.;
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      scalar hk = hl[nl-1-k];
      if (hk[] > dry) {
	H[l*nl + k] -= 8.*s*h[]/hk[];
	H[l*nl + k - 1] += 8.*s*h[]/hk[];
      }
    }
  }
}

#include "hessenberg.h"

trace
static void relax_nh (scalar * phil, scalar * rhsl, int lev, void * data)
{
  foreach_level_or_leaf (lev) {
    double H[nl*nl], b[nl];
    box_matrix (point, phil, rhsl, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
    for (scalar phi in phil)
      phi[] = b[l--];
  }
}

trace
static double residual_nh (scalar * phil, scalar * rhsl,
			   scalar * resl, void * data)
{
  double maxres = 0.;
#if 0 // TREE
  /* conservative coarse/fine discretisation (2nd order) */
  assert (false);
  face vector g[];
  foreach_face()
    g.x[] = face_gradient_x (a, 0);
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    scalar phi, rhs, res, h;
    int l = 0;
    for (phi,rhs,res,h in phil,rhsl,resl,hl) {
      res[] = rhs[] + 4.*phi[];
      foreach_dimension()
	res[] -= h[]*(h[1]*phi[1] - 2.*h[]*phi[] + h[-1]*phi[-1])/sq(Delta);
      if (l < nl - 1) {
        scalar phip = phil[l+1];
        res[] -= 4.*phip[];
	foreach_dimension()
	  res[] -= h[]*(h[1]*phip[1] - 2.*h[]*phip[] + h[-1]*phip[-1])/sq(Delta);
      }
      for (int k = l - 1, s = -1; k >= 0; k--, s = -s) {
	scalar hk = hl[k];
	if (hk[] > dry) {
	  scalar phik = phil[k], phikp = phil[k+1];
	  res[] += 8.*s*(phik[] - phikp[])*h[]/hk[];
	}
      }
#endif // !TREE    
#if EMBED
      if (p->embed_flux) {
	double c, e = p->embed_flux (point, a, alpha, &c);
	res[] += c - e*a[];
      }
#endif // EMBED    
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      l++;
    }
  }
  boundary (resl);
  return maxres;
}

double breaking = HUGE;

event pressure (i++)
{
  scalar * rhsl = list_clone (phil);
  foreach() {
    scalar rhs, h, qz;
    vector uf;
    int l = 0;
    coord zl, zr;
    foreach_dimension()
      zl.x = zb[-1], zr.x = zb[1];
    for (rhs,h,qz,uf in rhsl,hl,qzl,ufl) {
      rhs[] = 2.*qz[];
      foreach_dimension()
	rhs[] -= h[]*(uf.x[] + uf.x[1])*(zr.x + h[1] - zl.x - h[-1])/(4.*Delta);
      if (l > 0) {
	vector ufm = ufl[l-1];
	foreach_dimension()
	  rhs[] += h[]*(ufm.x[] + ufm.x[1])*(zr.x - zl.x)/(4.*Delta);
      }
      foreach_dimension()
	rhs[] += h[]*((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      for (int k = l - 1, s = -1; k >= 0; k--, s = -s) {
	scalar qk = qzl[k], hk = hl[k];
	if (hk[] > dry)
	  rhs[] += 4.*s*h[]*qk[]/hk[];
      }
      rhs[] *= 2./dt;
      foreach_dimension()
	zl.x += h[-1], zr.x += h[1];
      l++;
    }
  }
  
  restriction (hl);
  mgp = mg_solve (phil, rhsl, residual_nh, relax_nh, NULL,
		  nrelax = 4, res = NULL, minlevel = 1, tolerance = TOLERANCE);
  delete (rhsl), free (rhsl);
  
  foreach_face() {
    double H = 0., Hm = 0.;
    for (scalar h in hl)
      H += h[], Hm += h[-1];
    if ((H > dry && Hm > dry) ||
	(H > dry && eta[] >= zb[-1]) ||
	(Hm > dry && eta[-1] >= zb[])) {
      scalar phi, h;
      vector uf, a;
      int l = 0;
      double zr = zb[], zl = zb[-1];
      for (phi,h,uf,a in phil,hl,ufl,al) {
	if (h[] + h[-1] > dry) {
	  double ax;
#if 0
	  if (l == nl - 1)
	    ax = sq(fm.x[])*(h[]*phi[] - h[-1]*phi[-1] +
			     (phi[] + phi[-1])*(zr - zl))
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  else {
	    scalar phip = phil[l+1];
	    ax = sq(fm.x[])*(h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1])
			     - ((phip[] + phip[-1])*(zr + h[] - zl - h[-1]) -
				(phi[] + phi[-1])*(zr - zl)))
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  }
#else
	  if (l == nl - 1)
	    ax = fm.x[]*(h[]*phi[] - h[-1]*phi[-1] +
			 (phi[] + phi[-1])*(zr - zl))/(Delta*(h[] + h[-1]));
	  else {
	    scalar phip = phil[l+1];
	    ax = fm.x[]*(h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1])
			     - ((phip[] + phip[-1])*(zr + h[] - zl - h[-1]) -
				(phi[] + phi[-1])*(zr - zl)))
	      /(Delta*(h[] + h[-1]));
	  }
#endif
	  uf.x[] -= dt*ax;
	  a.x[] -= ax;
	}
	l++, zr += h[], zl += h[-1];
      }
    }
  }
  boundary ((scalar *)ufl);
  boundary_flux (al);

  foreach() {
    double wmax = 0.;
    for (scalar h in hl)
      wmax += h[];
    wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    scalar phi, qz, h;
    int l = 0;
    for (phi,qz,h in phil,qzl,hl) {
      if (l == nl - 1)
	qz[] += dt*phi[];
      else {
	scalar phip = phil[l+1];
	qz[] -= dt*(phip[] - phi[]);
      }
      if (fabs(qz[]) > h[]*wmax)
	qz[] = (qz[] > 0. ? 1. : -1.)*h[]*wmax;
      l++;
    }
  }
  boundary (qzl);
}

event cleanup (i = end, last) {
  delete (qzl), free (qzl), qzl = NULL;
  delete (phil), free (phil), phil = NULL;
}
