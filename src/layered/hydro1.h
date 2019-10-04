#include "run.h"
#if 1
#include "bcg.h"
#else
#include "mpdata.h"
#endif

scalar zb[], eta;
scalar * hl = NULL;
vector * ul = NULL;
double G = 1., dry = 1e-6;

int nl = 1;
double * beta = NULL;
double (* gradient) (double, double, double) = minmod2;

scalar ** tracers = NULL;

vector * ufl = NULL, * al = NULL;

attribute {
  int l;
}

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child() {
    eta[] = zb[];
    for (scalar h in hl)
      eta[] += h[];
  }
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[];
  for (scalar h in hl)
    eta[] += h[];
}

static const double default_sea_level = 0.;

static void refine_elevation (Point point, scalar h)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  if (h[] >= dry) {
    double eta = zb[] + h[];   // water surface elevation  
    coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
	g.x = gradient (zb[-1] + h[-1], eta, zb[1] + h[1])/4.;
    else
      foreach_dimension()
	g.x = (zb[1] - zb[-1])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta;
      foreach_dimension()
	etac += g.x*child.x;
      h[] = max(0., etac - zb[]);
    }
  }
  else {

    /**
    The "dry" case is a bit more complicated. We look in a 3x3
    neighborhood of the coarse parent cell and compute a depth-weighted
    average of the "wet" surface elevation $\eta$. We need to do this
    because we cannot assume a priori that the surrounding wet cells are
    necessarily close to e.g. $\eta = 0$. */

    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood
    foreach_neighbor(1)
      if (h[] >= dry) {
	eta += h[]*(zb[] + h[]);
	v += h[];
      }
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else

      /**
      If none of the surrounding cells is wet, we set a default sealevel. */

      eta = default_sea_level;

    /**
    We then reconstruct the water depth in each child using $\eta$ (of the
    parent cell i.e. a first-order interpolation in contrast to the wet
    case above) and $z_b$ of the child cells. */
    
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0., eta - zb[]);
  }
}

/**
Cell restriction is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */

static void restriction_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }

  /**
  ... and use this in combination with $z_b$ (of the coarse cell) to
  compute the water depth $h$.  */
    
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

/**
We also need to define a consistent prolongation function. For cells
which are entirely surrounded by wet cells, we can use the standard
linear refinement function, otherwise we use straight injection from
the parent cell. */

static void prolongation_elevation (Point point, scalar h)
{
  bool wet = true;
  foreach_neighbor(1)
    if (h[] <= dry)
      wet = false, break;
  if (wet)
    refine_linear (point, h);
  else {
    double hc = h[], zc = zb[];
    foreach_child() {
      h[] = hc;
      zb[] = zc;
    }
  }
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  for (scalar h in hl) {
    h.refine  = refine_elevation;
    h.prolongation = prolongation_elevation;
    h.restriction = restriction_elevation;
  }
}
#else
void conserve_elevation (void){}
#endif

event defaults0 (i = 0)
{
  assert (hl == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar h = new scalar;
    hl = list_append (hl, h);
    h.l = l;
  }
  eta = new scalar;
  reset (hl, 0.);
  reset ({zb}, 0.);

  if (!beta) {
    beta = malloc (nl*sizeof(double));
    for (int l = 0; l < nl; l++)
      beta[l] = 1./nl;
  }

  assert (!tracers);
  tracers = calloc (nl, sizeof(scalar *));

  int l = 0;
  for (scalar h in hl) {
    tracers[l] = list_prepend (tracers[l], h);
    l++;
  }

  zb.gradient = gradient;
#if TREE
  zb.refine = zb.prolongation = refine_linear;
  zb.restriction = restriction_volume_average;
  eta.prolongation = refine_linear;
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
  eta.gradient = gradient;
#endif
}

event defaults (i = 0)
{
  assert (ul == NULL && ufl == NULL && al == NULL);
  for (int l = 0; l < nl; l++) {
    vector u = new vector;
    vector uf = new face vector;
    vector a = new face vector;
    ul = vectors_append (ul, u);
    ufl = vectors_append (ufl, uf);
    al = vectors_append (al, a);
    foreach_dimension()
      u.x.l = uf.x.l = a.x.l = l;
  }
  reset (ul, 0.);
  reset (ufl, 0.);
  reset (al, 0.);

  int l = 0;
  for (vector u in ul) {
    foreach_dimension()
      tracers[l] = list_append (tracers[l], u.x);
    l++;
  }

  for (int l = 0; l < nl; l++)
    for (scalar s in tracers[l]) {
      s.gradient = gradient;
#if TREE
      s.refine = s.prolongation = refine_linear;
      s.restriction = restriction_volume_average;
#endif
    }
}

double dtmax;

event init (i = 0)
{
  trash ((scalar *)ufl);
  foreach_face() {
    vector uf, u;
    scalar h;
    for (h,u,uf in hl,ul,ufl)
#if 0      
      uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
#else
      uf.x[] = fm.x[]*(h[]*u.x[] + h[-1]*u.x[-1])/(h[] + h[-1] + dry);
#endif
  }
  boundary ((scalar *)ufl);
  dtmax = DT;
  event ("stability");
 
  foreach() {
    eta[] = zb[];
    for (scalar h in hl)
      eta[] += h[];
  }
  boundary (all);
}

event set_dtmax (i++,last) dtmax = DT;

static bool non_hydro = false;

scalar dtmin[];

event stability (i++,last)
{
#if TREE
  foreach_face() {
    double H = 0., Hm = 0.;
    scalar h;
    vector uf;
    for (h,uf in hl,ufl) {
      H += h[], Hm += h[-1];
      if ((h[] < dry && uf.x[] < 0.) ||
	  (h[-1] < dry && uf.x[] > 0.))
	uf.x[] = 0.;
    }
    if (!((H > dry && Hm > dry) ||
	  (H > dry && zb[] + H >= zb[-1]) ||
	  (Hm > dry && zb[-1] + Hm >= zb[]))) {
      for (vector uf in ufl)
	uf.x[] = 0.;
    }
  }
  boundary ((scalar *)ufl);
#endif

  reset ({dtmin}, nodata);
  Point q = {0};
  foreach_face(reduction(min:dtmax)) {
    double H = 0.;
    for (scalar h in hl)
      H += h[] + h[-1];
    if (H > 0.) {
      H /= 2.;
      double cp = non_hydro ? sqrt(G*Delta*tanh(H/Delta)) : sqrt(G*H);
      for (vector uf in ufl) {
	double c = fm.x[]*cp + fabs(uf.x[]);
	if (c > 0.) {
	  double dt = cm[]*Delta/c;
	  if (dt < dtmin[])
	    dtmin[] = dt;
	  if (dt < dtmax) {
	    dtmax = dt;
	    q = point;
	  }
	}
      }
    }
  }
  dt = dtnext (CFL*dtmax);
#if 0
  fprintf (stderr, "dt: %g %g %g\n", dt, dtmax, statsf(dtmin).min);
  debug (q);
#endif
}

static void advection1 (scalar * list, face vector uf, double dt)
{
  scalar h = list[0];
  assert (h.i == hl[h.l].i);
  face vector hflux[];
  tracer_fluxes (list[0], uf, hflux, dt, zeroc);
  list++;
  for (scalar f in list) {
    face vector flux[];
#if 1    
    foreach_face() {
      double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
      int i = -(s + 1.)/2.;
      double g = f.gradient ?
	f.gradient (f[i-1], f[i], f[i+1])/Delta :
	(f[i+1] - f[i-1])/(2.*Delta);
      double f2 = f[i] + s*(1. - s*un)*g*Delta/2.;
#if 0
      #if dimension > 1
      if (fm.y[i] && fm.y[i,1]) {
	double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
	double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
	f2 -= dt*vn*fyy/(2.*Delta);
      }
      #endif
      #if dimension > 2
      if (fm.z[i] && fm.z[i,0,1]) {
	double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
	double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
	f2 -= dt*wn*fzz/(2.*Delta);
      }
      #endif
#endif
      flux.x[] = f2*hflux.x[];
    }
    boundary_flux ({flux});
#else
    tracer_fluxes (f, uf, flux, dt, zeroc);
#endif
    foreach() {
      f[] *= h[];
      foreach_dimension()
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
    }
  }
  foreach() {
    foreach_dimension()
      h[] += dt*(hflux.x[] - hflux.x[1])/(Delta*cm[]);
    if (h[] < dry) {
      for (scalar f in list)
	f[] = 0.;
    }
    else
      for (scalar f in list)
	f[] /= h[];
  }
  boundary ({h});
  boundary (list);
}

event advection_term (i++,last)
{
  int l = 0;
  for (vector uf in ufl)
    advection1 (tracers[l++], uf, dt); // , (scalar *){g});

  foreach() {
    double H = 0.;
    for (scalar h in hl)
      H += h[];
    eta[] = zb[] + H;
  }
  boundary ({eta});
}

#include "viscosity-multilayer.h"

double filter = 0.;

scalar filtered[];

event viscous_term (i++,last)
{
#if 1
  if (filter > 0.) {
    assert (nl == 1); // fixme
    scalar h = hl[0];
    reset ({filtered}, 0);
    foreach()
      foreach_dimension()
        if (h[-1] > dry && h[] > dry && h[1] > dry) {
	  double Dp = eta[1] - eta[], Dm = eta[] - eta[-1];
	  if (Dp*Dm < 0. && ((eta[2] + eta[] - 2.*eta[1])*
			     (eta[-1] + eta[1] - 2.*eta[]) < 0. ||
			     (eta[-1] + eta[1] - 2.*eta[])*
			     (eta[-2] + eta[] - 2.*eta[-1]) < 0.)) {
	    double dp, dm;
	    if (fabs(Dp) > fabs(Dm)) {
	      dp = fabs(Dp);
	      dm = fabs(Dm);
	    }
	    else {
	      dp = fabs(Dm);
	      dm = fabs(Dp);
	    }
	    double d = min(dm, dp/2.);
	    double a = Dp > 0. ? 1. : -1.;
	    filtered[] = d;
#if 1	    
	    eta[] += min(dt/filter, 1.)*a*d;
	    h[] = max(eta[] - zb[], 0.);
            if (h[] < dry) // fixme for multiple layers
	      for (int l = 0; l < nl; l++) // fixme for multiple layers
		for (scalar s in tracers[l])
		  s[] = 0.;
#endif
	  }
	}
    boundary ({eta, h});
    for (int l = 0; l < nl; l++)
      boundary (tracers[l]);
  }
#endif  
  if (nu > 0.) {
#if 1
    foreach() {
#if 1 // does not seem very useful
      vector a, u;
      for (a,u in al,ul)
	foreach_dimension()
	  u.x[] += dt*(a.x[] + a.x[1])/(fm.x[] + fm.x[1] + SEPS);
#endif
      vertical_viscosity (point, hl, (scalar *) ul, dt); // fixme: 1D only
#if 1
      for (a,u in al,ul)
	foreach_dimension()
	  u.x[] -= dt*(a.x[] + a.x[1])/(fm.x[] + fm.x[1] + SEPS);
#endif
    }
#else
    foreach() {
      vector u;
      scalar h;
      int l = 0;
      for (u,h in ul,hl) {
	scalar up = l < nl - 1 ? ul[l+1].x : u.x;
	scalar hp = l < nl - 1 ? hl[l+1] : h;
	if (l == 0)
	  u.x[] += dt*nu*2.*((up[]/hp[] - u.x[]/h[])/(h[] + hp[]) -
			     u.x[]/su(h[]));
	else {
	  scalar um = ul[l-1].x, hm = hl[l-1];
	  u.x[] += dt*nu*2.*((up[]/hp[] - u.x[]/h[])/(h[] + hp[]) -
			     (u.x[]/h[] - um[]/hm[])/(h[] + hm[]));
	}
	l++;
      }
    }
#endif
    boundary ((scalar *) ul);

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
    vector u;
    if ((H > dry && Hm > dry) ||
	(H > dry && eta[] >= zb[-1]) ||
	(Hm > dry && eta[-1] >= zb[])) {
      for (h,a,uf,u in hl,al,ufl,ul) {
        a.x[] = - sq(fm.x[])*G*(eta[] - eta[-1])/((cm[] + cm[-1])*Delta/2.);
#if 0
        uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2. + dt*a.x[];
#else
        uf.x[] = fm.x[]*(h[]*u.x[] + h[-1]*u.x[-1])/(h[] + h[-1] + dry) +
	  dt*a.x[];
#endif
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
    vector a, u;
    for (a,u in al,ul) {
      foreach_dimension()
        u.x[] += dt*(a.x[] + a.x[1])/(fm.x[] + fm.x[1]);
#if dimension == 2
      // metric terms
      double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
      double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
      double ux = u.x[], uy = u.y[];
      double fG = uy*dmdl - ux*dmdt;
      u.x[] += dt*fG*uy;
      u.y[] -= dt*fG*ux;
#endif // dimension == 2
    }
  }
  boundary ((scalar *) ul);
}

#if TREE
event adapt (i++,last) {
}
#endif

event cleanup (i = end, last)
{
  free (beta), beta = NULL;
  for (int l = 0; l < nl; l++)
    free (tracers[l]);
  free (tracers), tracers = NULL;
  delete ((scalar *) ul), free (ul), ul = NULL;
  delete ((scalar *) ufl), free (ufl), ufl = NULL;
  delete ((scalar *) al), free (al), al = NULL;
  delete (hl), free (hl), hl = NULL;
  delete ({eta});
}

/**
# "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (*ref*). */

double _radiation (Point point, double ref, scalar s)
{
  double H = 0.;
  for (scalar h in hl)
    H += h[];
  return H > dry ? sqrt(G/H)*(zb[] + H - ref) : 0.;
}
  
#define radiation(ref) _radiation(point, ref, _s)

/**
# Tide gauges

An array of *Gauge* structures passed to *output_gauges()* will create
a file (called *name*) for each gauge. Each time *output_gauges()* is
called a line will be appended to the file. The line contains the time
and the value of each scalar in *list* in the cell containing
*(x,y)*. The *desc* field can be filled with a longer description of
the gauge. */

typedef struct {
  char * name;
  double x, y;
  char * desc;
  FILE * fp;
} Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  int n = 0;
  for (Gauge * g = gauges; g->name; g++, n++);
  coord a[n];
  n = 0;
  for (Gauge * g = gauges; g->name; g++, n++) {
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    a[n].x = xp, a[n].y = yp;
  }
  int len = list_len(list);
  double v[n*len];
  interpolate_array (list, a, n, v, false);

  if (pid() == 0) {
    n = 0;
    for (Gauge * g = gauges; g->name; g++) {
      if (!g->fp) {
	g->fp = fopen (g->name, "w");
	if (g->desc)
	  fprintf (g->fp, "%s\n", g->desc);
      }
      if (v[n] != nodata) {
	fprintf (g->fp, "%g", t);
	for (scalar s in list)
	  fprintf (g->fp, " %g", v[n++]);
	fputc ('\n', g->fp);
	fflush (g->fp);
      }
      else
	n += len;
    }
  }
}

void vertical_velocity (scalar * wl)
{
  foreach() {
    scalar h, w;
    vector uf;
    int l = 0;
    double zl = zb[-1], zr = zb[1];
    double wm = 0.;
    for (h,w,uf in hl,wl,ufl) {
      w[] = wm + (uf.x[] + uf.x[1])*(zr + h[1] - zl - h[-1])/(4.*Delta);
      if (l > 0) {
	vector ufm = ufl[l-1];
	w[] -= (ufm.x[] + ufm.x[1])*(zr - zl)/(4.*Delta);
      }
      foreach_dimension()
	w[] -= ((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      l++, zl += h[-1], zr += h[1], wm = w[];
    }
  }
  boundary (wl);
}
