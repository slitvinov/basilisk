A solver for the Saint-Venant equations
=======================================

First we need some utilities. Time integration will be done with a generic 
predictor-corrector scheme.

~~~
#include "utils.h"
#include "predictor-corrector.h"
~~~

The primary fields are the water depth h, the bathymetry zb and the flow 
speed u. eta is the water level i.e. h + zb. Note that the order of the 
declarations is important as zb needs to be refined before h and h before eta.

~~~
scalar zb[], h[], eta[];
vector u[];
~~~

The only physical parameter is the acceleration of gravity G. Cells are 
considered "dry" when the water depth is less than the dry parameter (this 
should not require tweaking). The initial conditions are defined by the 
user in init().

~~~
double G = 1.;
double dry = 1e-10;
void init (void);
~~~

The generic time-integration scheme in predictor-corrector.h needs to know 
which fields are updated.

~~~
scalar * evolving = {h, u};
~~~

When using an adaptive discretisation (i.e. a quadtree)., we need to make sure 
that eta is maintained as zb + h whenever cells are refined or coarsened.

~~~
#if QUADTREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h[];
}

static void coarsen_eta (Point point, scalar eta)
{
  eta[] = zb[] + h[];
}
#endif
~~~

The predictor-corrector scheme will call this function before starting the 
main time loop.

~~~
void init_internal (void)
{
#if QUADTREE
  eta.refine  = refine_eta;
  eta.coarsen = coarsen_eta;
#endif
  init();
  foreach()
    eta[] = zb[] + h[];
  boundary (all);
}
~~~

Computing fluxes
----------------

~~~
vector Fh, S;
tensor Fq;

#include "riemann.h"

double fluxes (scalar * evolving, double dtmax)
{
  // recover scalar and vector fields from list
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };

  // gradients
  vector gh[], geta[];
  tensor gu[];
  for (scalar s in {gh, geta, gu})
    s.gradient = zero; // first-order gradient reconstruction
  gradients ({h, eta, u}, {gh, geta, gu});

  // fluxes
  Fh = new vector; S = new vector;
  Fq = new tensor;
  foreach_face (reduction (min:dtmax)) {
    double hi = h[], hn = h[-1,0];
    if (hi > dry || hn > dry) {
      double dx = delta/2.;
      double zi = eta[] - hi;
      double zl = zi - dx*(geta.x[] - gh.x[]);
      double zn = eta[-1,0] - hn;
      double zr = zn + dx*(geta.x[-1,0] - gh.x[-1,0]);
      double zlr = max(zl, zr);
      
      double hl = hi - dx*gh.x[];
      double up = u.x[] - dx*gu.x.x[];
      double hp = max(0., hl + zl - zlr);
      
      double hr = hn + dx*gh.x[-1,0];
      double um = u.x[-1,0] + dx*gu.x.x[-1,0];
      double hm = max(0., hr + zr - zlr);

      // Riemann solver
      double fh, fu, fv;
      kurganov (hm, hp, um, up, delta, &fh, &fu, &dtmax);
      fv = (fh > 0. ? u.y[-1,0] + dx*gu.y.x[-1,0] : u.y[] - dx*gu.y.x[])*fh;

      // topographic source term
#if QUADTREE
      // this is necessary to ensure balance at fine/coarse boundaries
      // see notes/balanced.tm
      if (!(cell.flags & (fghost|active))) {
	hi = coarse(h,0,0);
	zi = coarse(zb,0,0);
      }
      if (!(neighbor(-1,0).flags & (fghost|active))) {
	hn = coarse(h,-1,0);
	zn = coarse(zb,-1,0);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

      // fluxes
      Fh.x[]   = fh;
      Fq.x.x[] = fu - sl;
      S.x[]    = fu - sr;
      Fq.y.x[] = fv;
    }
    else // dry
      Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
  }

#if QUADTREE
  // fixme: all this should be halo_restriction_flux()
  vector * list = {Fh, S, Fq};
  foreach_halo_fine_to_coarse()
    foreach_dimension() {
      if (is_leaf (neighbor(-1,0)))
	for (vector f in list)
	  f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
      if (is_leaf (neighbor(1,0)))
	for (vector f in list)
	  f.x[1,0] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
    }
#endif

  return dtmax;
}

void update (scalar * output, scalar * input, double dt)
{
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0];
  vector 
    ui = { input[1], input[2] },
    uo = { output[1], output[2] };

  if (ho != hi)
    trash ({ho, uo});
  trash ({eta});
  // new fields in ho[], uo[]
  foreach() {
    double hold = hi[];
    ho[] = hold + dt*(Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/delta;
    eta[] = ho[] + zb[];
    if (ho[] > dry)
      foreach_dimension()
	uo.x[] = (hold*ui.x[] + dt*(Fq.x.x[] + Fq.x.y[] - 
				   S.x[1,0] - Fq.x.y[0,1])/delta)/ho[];
    else
      foreach_dimension()
	uo.x[] = 0.;
  }
  boundary ({ho, eta, uo});

  delete ((scalar *){Fh, S, Fq});
}

#if QUADTREE
static void refine_elevation (Point point, scalar h)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  if (h[] >= dry) {
    // wet cell
    double eta = zb[] + h[];   // water surface elevation  
    struct { double x, y; } g; // gradient of eta
    foreach_dimension()
      g.x = gradient (zb[-1,0] + h[-1,0], eta, zb[1,0] + h[1,0])/4.;
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta + g.x*child.x + g.y*child.y - zb[]);
  }
  else {
    // dry cell
    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h[i,j] >= dry) {
	  eta += h[i,j]*(zb[i,j] + h[i,j]);
	  v += h[i,j];
	}
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else
      eta = 0.; // surrounded by dry cells => set eta to default "sealevel"

    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta - zb[]);
  }
}

static void coarsen_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

void conserve_elevation (void)
{
  h.refine  = refine_elevation;
  h.coarsen = coarsen_elevation;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif

// "radiation" boundary conditions
#define radiation(ref) (sqrt (G*h[]) - sqrt(G*max((ref) - zb[], 0.)))

// tide gauges

typedef struct {
  char * name;
  double x, y;
  char * desc;
  FILE * fp;
} Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  for (Gauge * g = gauges; g->name; g++) {
    if (!g->fp) {
      g->fp = fopen (g->name, "w");
      if (g->desc)
	fprintf (g->fp, "%s\n", g->desc);
    }
    Point point = locate (g->x, g->y);
    if (point.level >= 0 && h[] > dry) {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " %g", s[]);
    }
    fputc ('\n', g->fp);
  }
}
~~~
