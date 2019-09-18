#include "run.h"
#include "bcg.h"

scalar zb[];
scalar * hl = NULL;
vector * ql = NULL;
double G = 1., dry = 1e-6;

int nl = 1;
double * beta = NULL;

scalar ** tracers = NULL;

vector * ufl = NULL, * al = NULL;

attribute {
  int l;
}

event defaults (i = 0)
{
  assert (hl == NULL && ql == NULL && ufl == NULL && al == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    vector q = new vector;
    vector uf = new face vector;
    scalar h = new scalar;
    vector a = new face vector;
    ql = vectors_append (ql, q);
    ufl = vectors_append (ufl, uf);
    hl = list_append (hl, h);
    al = vectors_append (al, a);
    foreach_dimension()
      q.x.l = l;
  }
  reset ({zb}, 0.);
  reset (hl, 0.);
  reset (ql, 0.);
  reset (al, 0.);

  if (!beta) {
    beta = malloc (nl*sizeof(double));
    for (int l = 0; l < nl; l++)
      beta[l] = 1./nl;
  }

  if (!tracers)
    tracers = calloc (nl, sizeof(scalar *));

  scalar h;
  vector q;
  int l = 0;
  for (h,q in hl,ql) {
    tracers[l] = list_prepend (tracers[l], h);
    foreach_dimension()
      tracers[l] = list_append (tracers[l], q.x);
    l++;
  }
}

double dtmax;

event init (i = 0)
{
  boundary (all);
  
  trash ((scalar *)ufl);
  foreach_face() {
    vector uf, q;
    scalar h;
    for (h,q,uf in hl,ql,ufl)
      uf.x[] = fm.x[]*(q.x[] + q.x[-1])/(h[] + h[-1] + dry);
  }
  boundary ((scalar *)ufl);
  dtmax = DT;
  event ("stability");

  zb.gradient = minmod2;
  for (scalar s in hl)
    s.gradient = minmod2;
  for (vector q in ql)
    foreach_dimension()
      q.x.gradient = minmod2;
  //  theta = 1.;
}

event set_dtmax (i++,last) dtmax = DT;

static bool non_hydro = false;

event stability (i++,last)
{
  foreach_face(reduction(min:dtmax)) {
    double H = 0.;
    for (scalar h in hl)
      H += h[] + h[-1];
    H /= 2.;
    double cp = non_hydro ? sqrt(G*Delta*tanh(H/Delta)) : sqrt(G*H);
    for (vector uf in ufl) {
      double c = cp + fabs(uf.x[]);
      if (c > 0.) {
	double dt = Delta/c;
	if (dt < dtmax)
	  dtmax = dt;
      }
    }
  }
  dt = dtnext (CFL*dtmax);
}

event advection_term (i++,last)
{
  int l = 0;
  for (vector uf in ufl)
    advection (tracers[l++], uf, dt); // , (scalar *){g});
}

#include "viscosity-multilayer.h"

event viscous_term (i++,last)
{
  if (nu > 0.) {
#if 1
    foreach() {
#if 1 // does not seem very useful
      scalar h;
      vector a, q;
      for (a,h,q in al,hl,ql)
	foreach_dimension()
	  q.x[] += dt*h[]*(a.x[] + a.x[1])/(fm.x[] + fm.x[1] + SEPS);
#endif
      vertical_viscosity (point, hl, (scalar *) ql, dt); // fixme: 1D only
#if 1
      for (a,h,q in al,hl,ql)
	foreach_dimension()
	  q.x[] -= dt*h[]*(a.x[] + a.x[1])/(fm.x[] + fm.x[1] + SEPS);
#endif
    }
#else
    foreach() {
      vector q;
      scalar h;
      int l = 0;
      for (q,h in ql,hl) {
	scalar qp = l < nl - 1 ? ql[l+1].x : q.x;
	scalar hp = l < nl - 1 ? hl[l+1] : h;
	if (l == 0)
	  q.x[] += dt*nu*2.*((qp[]/hp[] - q.x[]/h[])/(h[] + hp[]) -
			     q.x[]/sq(h[]));
	else {
	  scalar qm = ql[l-1].x, hm = hl[l-1];
	  q.x[] += dt*nu*2.*((qp[]/hp[] - q.x[]/h[])/(h[] + hp[]) -
			     (q.x[]/h[] - qm[]/hm[])/(h[] + hm[]));
	}
	l++;
      }
    }
#endif
    boundary ((scalar *) ql);
  }
}

event acceleration (i++,last)
{
  trash (ufl);
  foreach_face() {
    double H = 0., Hm = 0.;
    for (scalar h in hl)
      H += h[], Hm += h[-1];
    scalar h;
    face vector a, uf;
    vector q;
    if ((H > dry && Hm > dry) ||
	(H > dry && H + zb[] >= zb[-1]) ||
	(Hm > dry && Hm + zb[-1] >= zb[])) {
      for (h,a,uf,q in hl,al,ufl,ql) {
	if (h[] > dry && h[-1] > dry) {
	  a.x[] = - fm.x[]*G*(zb[] + H - zb[-1] - Hm)/Delta;
	  uf.x[] = fm.x[]*(q.x[] + q.x[-1])/(h[] + h[-1]) + dt*a.x[];
	}
	else
	  a.x[] = uf.x[] = 0.;
      }
    }
    else {
      for (a,uf in al,ufl)
	a.x[] = uf.x[] = 0.;
    }
  }
  boundary ((scalar *)ufl);
  boundary ((scalar *)al);
}

event pressure (i++,last)
{
  foreach() {
    scalar h;
    vector a, q;
    for (a,h,q in al,hl,ql)
      foreach_dimension()
	q.x[] += dt*h[]*(a.x[] + a.x[1])/(fm.x[] + fm.x[1] + SEPS);
  }
  boundary ((scalar *) ql);
}

#if TREE
event adapt (i++,last);
#endif

event cleanup (i = end, last) {
  free (beta), beta = NULL;
  for (int l = 0; l < nl; l++)
    free (tracers[l]);
  free (tracers), tracers = NULL;
  delete ((scalar *) ql), free (ql), ql = NULL;
  delete ((scalar *) ufl), free (ufl), ufl = NULL;
  delete ((scalar *) al), free (al), al = NULL;
  delete (hl), free (hl), hl = NULL;
}
