#define NH 1
#define BOX 1 // fixme

#include "poisson.h"

scalar * wl = NULL, * phil = NULL;
mgstats mgp;

event defaults (i = 0)
{
  non_hydro = true;
  mgp.nrelax = 4;
  
  assert (wl == NULL && phil == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar w = new scalar;
    scalar phi = new scalar;
    wl = list_append (wl, w);
    phil = list_append (phil, phi);
  }
  reset (wl, 0.);
  reset (phil, 0.);

#if 0  
  if (!tracers)
    tracers = calloc (nl, sizeof(scalar *));
#endif
  
  int l = 0;
  for (scalar w in wl) {
    tracers[l] = list_append (tracers[l], w);
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
      vertical_viscosity (point, hl, (scalar *) wl, dt);
    boundary (wl);
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

#define slope_limited(dz) (fabs(dz) < 1. ? (dz) : ((dz) > 0. ? 1. : -1.))
 
event pressure (i++)
{
  double h1 = 0., v1 = 0.;
  scalar * rhsl = list_clone (phil);
  foreach (reduction(+:h1) reduction(+:v1)) {
    scalar rhs, h, w;
    vector uf;
    int l = 0;
    coord dz;
    foreach_dimension()
      dz.x = zb[1] - zb[-1];
    for (rhs,h,w,uf in rhsl,hl,wl,ufl) {
      rhs[] = 2.*h[]*w[];
      foreach_dimension()
	rhs[] -= h[]*(uf.x[] + uf.x[1])/2.*
	slope_limited((dz.x + h[1] - h[-1])/(2.*Delta));
      if (l > 0) {
	vector ufm = ufl[l-1];
	foreach_dimension()
	  rhs[] += h[]*(ufm.x[] + ufm.x[1])/2.*slope_limited(dz.x/(2.*Delta));
      }
      foreach_dimension()
	rhs[] += h[]*((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      for (int k = l - 1, s = -1; k >= 0; k--, s = -s) {
	scalar uk = wl[k];
	rhs[] += 4.*s*h[]*uk[];
      }
      rhs[] *= 2./dt;
      foreach_dimension()
	dz.x += h[1] - h[-1];
      h1 += dv()*h[];
      v1 += dv();
      l++;
    }
  }
  
  restriction (hl);
  mgp = mg_solve (phil, rhsl, residual_nh, relax_nh, NULL,
		  nrelax = 4, res = NULL, minlevel = 1,
#if 0
		  tolerance = TOLERANCE
#else
		  tolerance = TOLERANCE*sq(h1/(dt*v1))
#endif
		  );
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
      double dz = zb[] - zb[-1];
      for (phi,h,uf,a in phil,hl,ufl,al) {
	if (h[] + h[-1] > dry) {
	  double ax;
#if 0
	  if (l == nl - 1)
	    ax = sq(fm.x[])*(h[]*phi[] - h[-1]*phi[-1] +
			     (phi[] + phi[-1])*dz)
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  else {
	    scalar phip = phil[l+1];
	    ax = sq(fm.x[])*(h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1])
			     - ((phip[] + phip[-1])*(dz + h[] - h[-1]) -
				(phi[] + phi[-1])*dz))
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  }
#else
	  if (l == nl - 1)
	    ax = fm.x[]*((h[]*phi[] - h[-1]*phi[-1])/Delta +
			 (phi[] + phi[-1])*slope_limited(dz/Delta))
	      /(h[] + h[-1]);
	  else {
	    scalar phip = phil[l+1];
	    ax = fm.x[]*
	      ((h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1]))/Delta
	       - ((phip[] + phip[-1])*slope_limited((dz + h[] - h[-1])/Delta) -
		  (phi[] + phi[-1])*slope_limited(dz/Delta)))
	      /(h[] + h[-1]);
	  }
#endif
	  uf.x[] -= dt*ax;
	  a.x[] -= ax;
	}
	l++, dz += h[] - h[-1];
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
    scalar phi, w, h;
    int l = 0;
    for (phi,w,h in phil,wl,hl) {
      if (h[] > dry) {
	if (l == nl - 1)
	  w[] += dt*phi[]/h[];
	else {
	  scalar phip = phil[l+1];
	  w[] -= dt*(phip[] - phi[])/h[];
	}
	if (fabs(w[]) > wmax)
	  w[] = (w[] > 0. ? 1. : -1.)*wmax;
      }
      l++;
    }
  }
  boundary (wl);
}

event cleanup (i = end, last) {
  delete (wl), free (wl), wl = NULL;
  delete (phil), free (phil), phil = NULL;
}
