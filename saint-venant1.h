#include "utils.h"

// h: water depth
// zb: bathymetry
// u: flow speed
scalar h[], zb[];
vector u[];

// Default parameters
// acceleration of gravity
double G = 1.;
// dry
double dry = 1e-10;

// fields updated by time-integration
scalar * evolving = {h, u};

#include "riemann.h"

double tendencies (scalar * evolution, scalar * evolving, double dtmax)
{
  // recover scalar and vector fields from lists
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };
  scalar dh = evolution[0];
  vector dq = { evolution[1], evolution[2] };

  // gradients
  vector gh[], gzb[];
  tensor gu[];
  // first-order gradient reconstruction (for consistent limiting)
  for (scalar s in {gh, gzb, gu})
    s.gradient = zero;
  gradients ({h, zb, u}, {gh, gzb, gu});

  // fluxes
  foreach_face (reduction (min:dtmax)) {
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
      if (!(neighbor(-1,0).flags & (fghost|active))) {
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
  // propagate tendencies from fine to coarse
  foreach_halo()
    for (scalar ds in evolution) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = 0.;
    }
#endif

  return dtmax;
}

void update (scalar * evolving1, 
	     scalar * evolving, scalar * evolution, double dt)
{
  // recover scalar and vector fields from lists
  scalar h = evolving[0], h1 = evolving1[0];
  vector u = { evolving[1], evolving[2] }, u1 = { evolving1[1], evolving1[2] };
  scalar dh = evolution[0];
  vector dq = { evolution[1], evolution[2] };

  if (evolving1 != evolving)
    trash (evolving1);
  foreach() {
    double hold = h[];
    h1[] = hold + dt*dh[]/delta;
    dh[] = 0.;
    if (h1[] > dry)
      foreach_dimension() {
	u1.x[] = (hold*u.x[] + dt*dq.x[]/delta)/h1[];
	dq.x[] = 0.;
      }
    else
      foreach_dimension()
	u1.x[] = dq.x[] = 0.;
  }
  boundary (evolving1);
}

#include "predictor-corrector.h"
