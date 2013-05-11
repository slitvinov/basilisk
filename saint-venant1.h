#include "utils.h"

// h: water depth
// zb: bathymetry
// u: flow speed
scalar h = new scalar, zb = new scalar;
vector u = new vector;
// extra storage for predictor/corrector
// fixme: not needed for first-order scheme
scalar h1 = new scalar;
vector u1 = new vector;
// gradients
vector gh = new vector, gzb = new vector;
tensor gu = new tensor;
// tendencies
scalar dh = new scalar;
vector dq = new vector;

// Default parameters
// acceleration of gravity
double G = 1.;
// gradient
void (* gradient) (scalar *, vector *) = generalized_minmod;
// dry
double dry = 1e-10;
// user-provided functions
void parameters (void);
void init       (void);

#include "riemann.h"

static double flux (double dtmax)
{
  foreach_face() {
    double hi = h[], hn = h[-1,0];
    if (hi > dry || hn > dry) {
      double dx = delta/2.;
      double zl = zb[] - dx*gzb.x[];
      double zr = zb[-1,0] + dx*gzb.x[-1,0];
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
      double zi = zb[], zn = zb[-1,0];
#if QUADTREE
      // this is necessary to ensure balance at fine/coarse boundaries
      // see notes/balanced.tm
      if (!(cell.flags & (fghost|active))) {
	hi = coarse(h,0,0);
	zi = coarse(zb,0,0);
      }
      if (!(neighbor(-1,0).flags &(fghost|active))) {
	hn = coarse(h,-1,0);
	zn = coarse(zb,-1,0);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

      // update tendencies
        dh[] += fh;         dh[-1,0] -= fh;
      dq.y[] += fv;       dq.y[-1,0] -= fv;
      dq.x[] += fu - sl;  dq.x[-1,0] -= fu - sr;
    }
  }

#if QUADTREE
  // propagate updates from fine to coarse
  scalar * list = scalars(dh, dq);
  foreach_halo()
    for (scalar ds in list) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = 0.;
    }
#endif

  return dtmax;
}

static void update (vector u2, vector u1, scalar h2, scalar h1, double dt)
{
  if (h1 != h2)
    trash (h1, u1);
  foreach() {
    h1[] = h2[] + dt*dh[]/delta;
    dh[] = 0.;
    if (h1[] > dry)
      foreach_dimension() {
	u1.x[] = (h2[]*u2.x[] + dt*dq.x[]/delta)/h1[];
	dq.x[] = 0.;
      }
    else
      foreach_dimension()
	u1.x[] = dq.x[] = 0.;
  }
  boundary (h1, u1);
}

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

#if QUADTREE
  // we need the tendencies to be reinitialised during refinement
  scalar * tendencies = scalars (dh, dq);
  for (scalar ds in tendencies)
    _refine[ds] = refine_reset;
#endif

  // default values
  scalar * list = scalars (h, zb, u, dh, dq);
  foreach() {
    for (scalar s in list)
      s[] = 0.;
    h[] = 1.;
  }
  boundary (list);

  // user-defined initial conditions
  init();
  boundary (list);

  // clone temporary storage
  clone_scalar (h, h1);
  foreach_dimension()
    clone_scalar (u.x, u1.x);

  // main loop
  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    gradient (scalars (h, zb, u), vectors (gh, gzb, gu));
    dt = dtnext (t, flux (DT));

    if (gradient == zero)
      /* 1st-order time-integration for 1st-order spatial
	 discretisation */
      update (u, u, h, h, dt);
    else {
      /* 2nd-order time-integration */
      /* predictor */
      update (u, u1, h, h1, dt/2.);
      swap (vector, u, u1);
      swap (scalar, h, h1);
      
      /* corrector */
      gradient (scalars (h, u), vectors (gh, gu));
      flux (dt);

      update (u1, u, h1, h, dt);
    }

    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
