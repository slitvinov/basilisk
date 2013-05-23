#include "utils.h"
#include "predictor-corrector.h"

// h: water depth
// zb: bathymetry
// u: flow speed

 // note: the order is important as zb needs to be refined before h
scalar zb[], h[];
vector u[];

// Default parameters
// acceleration of gravity
double G = 1.;
// cells are "dry" when water depth is less than...
double dry = 1e-10;

// fields updated by time-integration (in "predictor-corrector.h")
scalar * evolving = {h, u};

#include "riemann.h"

// fluxes
vector Fh, S;
tensor Fq;

double fluxes (scalar * evolving, double dtmax)
{
  // recover scalar and vector fields from lists
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };

  // gradients
  vector gh[], gzb[];
  tensor gu[];
  for (scalar s in {gh, gzb, gu})
    s.gradient = zero;
  gradients ({h, u}, {gh, gu});
  // reconstruct zb + h rather than zb: see theorem 3.1 of Audusse et al, 2004
  foreach()
    foreach_dimension()
      gzb.x[] = gradient (zb[-1,0] + h[-1,0], 
			  zb[] + h[], 
			  zb[1,0] + h[1,0])/delta - gh.x[];
  boundary ((scalar *) {gzb});

  // fluxes
  Fh = new vector; S = new vector;
  Fq = new tensor;
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
  scalar h = input[0], ho = output[0];
  vector 
    u = { input[1], input[2] },
    uo = { output[1], output[2] };

  if (ho != h)
    trash ({ho, uo});
  // new fields in ho[], uo[]
  foreach() {
    double hold = h[];
    ho[] = hold + dt*(Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/delta;
    if (ho[] > dry)
      foreach_dimension()
	uo.x[] = (hold*u.x[] + dt*(Fq.x.x[] + Fq.x.y[] - 
				   S.x[1,0] - Fq.x.y[0,1])/delta)/ho[];
    else
      foreach_dimension()
	uo.x[] = 0.;
  }
  boundary (output);

  delete ((scalar *){Fh, S, Fq});
}
