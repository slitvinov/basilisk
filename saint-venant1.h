#include "utils.h"

// h: water depth
// zb: bathymetry
// q: flow rate
scalar h = new scalar, zb = new scalar;
vector q = new vector;
// storage for predictor/corrector
scalar h1 = new scalar;
vector q1 = new vector;
// topographic source term
vector Sb = new vector;
// gradients
vector gh = new vector, gzb = new vector;
tensor gq = new tensor;
// fluxes
vector fh = new vector;
tensor fq = new tensor;

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
    double eta = h[], etan = h[-1,0];
    if (eta <= dry && etan <= dry)
      fh.x[] = fq.x.x[] = fq.y.x[] = 0.;
    else {
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
      kurganov (hR, hL, uR, uL, DX, &fh.x[], &fq.x.x[], &dtmax);
      fq.y.x[] = (fh.x[] > 0. ? vR : vL)*fh.x[];

      // topographic source term
      if (eta <= dry) eta = 0.;
      if (etan <= dry) etan = 0.;
      Sb.x[] -= G/2.*(sq(hL) - sq(etaL) + dx*(etaL + eta)*gzb.x[]);
      Sb.x[-1,0] += G/2.*(sq(hR) - sq(etaR) - dx*(etaR + etan)*gzb.x[-1,0]);
    }
  }
  boundary_flux (fh, fq.x, fq.y);
  return dtmax;
}

static void update (vector q2, vector q1, scalar h2, scalar h1, double dt)
{
  foreach() {
    h1[] = h2[] + dt*(fh.x[] - fh.x[1,0] + fh.y[] - fh.y[0,1])/DX;
    foreach_dimension() {
      q1.x[] = q2.x[] + dt*(fq.x.x[] - fq.x.x[1,0] + fq.x.y[] - fq.x.y[0,1] 
			    + Sb.x[])/DX;
      Sb.x[] = 0.;
    }
  }
  boundary (h1, q1.x, q1.y);
}

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);
  init();
  boundary (h, zb, q.x, q.y);

  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    (* gradient) (scalars (q.x, q.y, h, zb), vectors (gq.x, gq.y, gh, gzb));
    boundary (gq.x.x, gq.x.y, gq.y.x, gq.y.y, gh.x, gh.y, gzb.x, gzb.y);

    dt = dtnext (t, flux (DT));

    if (gradient == zero)
      /* 1st-order time-integration for 1st-order spatial
	 discretisation */
      update (q, q, h, h, dt);
    else {
      /* 2nd-order time-integration */
      /* predictor */
      update (q, q1, h, h1, dt/2.);
      swap (vector, q, q1);
      swap (scalar, h, h1);
      
      /* corrector */
      (* gradient) (scalars (q.x, q.y, h), vectors (gq.x, gq.y, gh));
      boundary (gq.x.x, gq.x.y, gq.y.x, gq.y.y, gh.x, gh.y);
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
