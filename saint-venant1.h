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

double fluxes (scalar * evolving, vector Fh, vector S, tensor Fq, double dtmax)
{
  // recover scalar and vector fields from lists
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };

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

void update (scalar * evolving1, scalar * evolving, 
	     vector Fh, vector S, tensor Fq, double dt)
{
  // recover scalar and vector fields from lists
  scalar h = evolving[0], h1 = evolving1[0];
  vector 
    u = { evolving[1], evolving[2] },
    u1 = { evolving1[1], evolving1[2] };

  if (h1 != h)
    trash ({h1, u1});
  foreach() {
    double hold = h[];
    h1[] = hold + dt*(Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/delta;
    if (h1[] > dry)
      foreach_dimension()
	u1.x[] = (hold*u.x[] + dt*(Fq.x.x[] + Fq.x.y[] - 
				   S.x[1,0] - Fq.x.y[0,1])/delta)/h1[];
    else
      foreach_dimension()
	u1.x[] = 0.;
  }
  boundary (evolving1);
}

// Generic predictor/corrector time-integration

// User-provided parameters/functions
// gradient
double (* gradient)  (double, double, double) = minmod2;
void      parameters (void);
void      init       (void);

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

  // limiting
  for (scalar s in all)
    s.gradient = gradient;

  // default values
  foreach()
    for (scalar s in all)
      s[] = 0.;
  boundary (all);

  // user-defined initial conditions
  init();
  boundary (all);

  // main loop
  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    vector Fh[], S[];
    tensor Fq[];
    dt = dtnext (t, fluxes (evolving, Fh, S, Fq, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      // temporary storage
      scalar * temporary = clone (evolving);
      /* predictor */
      update (temporary, evolving, Fh, S, Fq, dt/2.);
      /* corrector */
      fluxes (temporary, Fh, S, Fq, dt);
      // free temporary storage
      delete (temporary);
      free (temporary);
    }
    update (evolving, evolving, Fh, S, Fq, dt);

    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
