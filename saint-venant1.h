#include "utils.h"

// h: water depth
// zb: bathymetry
// q: flow rate
scalar h = new scalar, zb = new scalar;
vector q = new vector;
// storage for predictor/corrector
scalar h1 = new scalar;
vector q1 = new vector;
// gradients
vector gh = new vector, gzb = new vector;
tensor gq = new tensor;
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

#if !QUADTREE
# define trash(x)
#endif

static double flux (double dtmax)
{
  trash (dh, dq);

  scalar * list = scalars(dh, dq);
#if QUADTREE
  foreach_halo()
    for (scalar ds in list) {
      coarse(ds,0,0) = 0.;
      ds[] = 0.;
    }
#endif
  foreach_face()
    for (scalar ds in list)
      ds[] = ds[-1,0] = 0.;

  foreach_face() {
    double eta = h[], etan = h[-1,0];
    if (eta > dry || etan > dry) {
      double dx = delta/2.;
      double zbL = zb[] - dx*gzb.x[];
      double zbR = zb[-1,0] + dx*gzb.x[-1,0];
      double zbLR = max(zbL, zbR);
      
      double etaL = eta <= dry ? 0. : eta - dx*gh.x[];
      double uL, vL;
      if (etaL > dry) {
	uL = (q.x[] - dx*gq.x.x[])/etaL;
	vL = (q.y[] - dx*gq.y.x[])/etaL;
      }
      else
	uL = vL = 0.;
      double hL = max(0., etaL + zbL - zbLR);
      
      double etaR = etan <= dry ? 0. : etan + dx*gh.x[-1,0];
      double uR, vR;
      if (etaR > dry) {
	uR = (q.x[-1,0] + dx*gq.x.x[-1,0])/etaR;
	vR = (q.y[-1,0] + dx*gq.y.x[-1,0])/etaR;
      }
      else
	uR = vR = 0.;
      double hR = max(0., etaR + zbR - zbLR);

      // Riemann solver
      double fh, fu, fv;
      kurganov (hR, hL, uR, uL, DX, &fh, &fu, &dtmax);
      fv = (fh > 0. ? vR : vL)*fh;

      // topographic source term
      if (eta <= dry) eta = 0.;
      if (etan <= dry) etan = 0.;
      double SbL = G/2.*(sq(hL) - sq(etaL) + dx*(etaL + eta)*gzb.x[]);
      double SbR = G/2.*(sq(hR) - sq(etaR) - dx*(etaR + etan)*gzb.x[-1,0]);

      // update tendencies
        dh[] += fh;           dh[-1,0] -= fh;
      dq.y[] += fv;         dq.y[-1,0] -= fv;
      dq.x[] += fu - SbL;   dq.x[-1,0] -= fu - SbR;
    }
  }

#if QUADTREE
  foreach_halo()
    for (scalar ds in list) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = undefined;
    }
#endif

  return dtmax;
}

static void update (vector q2, vector q1, scalar h2, scalar h1, double dt)
{
  trash (h1, q1);
  foreach() {
    h1[] = h2[] + dt*dh[]/DX;
    dh[] = undefined;
    foreach_dimension() {
      q1.x[] = q2.x[] + dt*dq.x[]/DX;
      dq.x[] = undefined;
    }
  }
  boundary (h1, q1);
}

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

  foreach() {
    zb[] = q.x[] = q.y[] = 0.;
    h[] = 1.;
  }
  init();
  boundary (h, zb, q);

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    (* gradient) (scalars (h, zb, q), vectors (gh, gzb, gq));

    dt = dtnext (t, flux (DT));

    if (gradient == zero)
      /* 1st-order time-integration for 1st-order spatial
	 discretisation */
      update (q, q, h, h, dt);
    else {
      /* 2nd-order time-integration */
      /* predictor */
      update (q, q1, h, h1, dt/2.);
      swap (vector, q, q1); // fixme: boundary conditions?
      swap (scalar, h, h1);
      
      /* corrector */
      (* gradient) (scalars (h, q), vectors (gh, gq));
      flux (dt);
      update (q1, q, h1, h, dt);
    }

    foreach (reduction(+:tnc)) 
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
