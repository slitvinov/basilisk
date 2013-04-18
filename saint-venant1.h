#include "utils.h"
#include "events.h"

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

#include "boundary.h"

// Default parameters
// acceleration of gravity
double G = 1.;
// gradient
void (* gradient) (const scalar, vector) = generalized_minmod;
// dry
double dry = 1e-10;
// user-provided functions
void parameters (void);
void init       (void);

#include "riemann.h"

static double flux (Point point, int i, double dtmax)
{
  foreach_dimension() {
    VARIABLES;
    double eta = h[i,0], etan = h[i-1,0];
    if (eta <= dry && etan <= dry)
      fh.x[i,0] = fq.x.x[i,0] = fq.y.x[i,0] = 0.;
    else {
      double zbL = zb[i,0] - gzb.x[i,0]/2.;
      double zbR = zb[i-1,0] + gzb.x[i-1,0]/2.;
      double zbLR = max(zbL, zbR);
      
      double etaL = eta <= dry ? 0. : eta - gh.x[i,0]/2.;
      double uL, vL;
      if (etaL > dry) {
	uL = (q.x[i,0] - gq.x.x[i,0]/2.)/etaL;
	vL = (q.y[i,0] - gq.y.x[i,0]/2.)/etaL;
      }
      else
	uL = vL = 0.;
      double hL = max(0., etaL + zbL - zbLR);
      
      double etaR = etan <= dry ? 0. : etan + gh.x[i-1,0]/2.;
      double uR, vR;
      if (etaR > dry) {
	uR = (q.x[i-1,0] + gq.x.x[i-1,0]/2.)/etaR;
	vR = (q.y[i-1,0] + gq.y.x[i-1,0]/2.)/etaR;
      }
      else
	uR = vR = 0.;
      double hR = max(0., etaR + zbR - zbLR);
      
      kurganov (hR, hL, uR, uL, DX, &fh.x[i,0], &fq.x.x[i,0], &dtmax);
      fq.y.x[i,0] = (fh.x[i,0] > 0. ? vR : vL)*fh.x[i,0];

      /* topographic source term */
      if (eta <= dry) eta = 0.;
      if (etan <= dry) etan = 0.;
      Sb.x[i,0] -= G/2.*(sq(hL) - sq(etaL) + (etaL + eta)*gzb.x[i,0]/2.);
      Sb.x[i-1,0] += G/2.*(sq(hR) - sq(etaR) - (etaR + etan)*gzb.x[i-1,0]/2.);
    }
  }
  return dtmax;
}

#define vswap(a,b) { vector tmp = a; a = b; b = tmp; }
#define swap(a,b) { scalar tmp = a; a = b; b = tmp; }

static void gradients ()
{
  (* gradient) (q.x, gq.x); boundary (gq.x.x); boundary (gq.x.y);
  (* gradient) (q.y, gq.y); boundary (gq.y.x); boundary (gq.y.y);
  (* gradient) (h, gh); boundary (gh.x); boundary (gh.y);
}

static double fluxes (double dt)
{
  foreach_boundary (right)
    dt = flux (point, 1, dt);
  foreach_boundary (top)
    dt = flux (point, 1, dt);
  foreach (reduction(min:dt))
    dt = flux (point, 0, dt);
  return dt;
}

static void update (vector q2, vector q1, scalar h2, scalar h1, double dt)
{
  foreach() {
    h1[] = h2[] + dt*(fh.x[] - fh.x[1,0] + fh.y[] - fh.y[0,1])/DX;
    foreach_dimension() {
      q1.x[] = q2.x[] + dt*(fq.x.x[] - fq.x.x[1,0] + fq.x.y[] - fq.x.y[0,1] + Sb.x[])/DX;
      Sb.x[] = 0.;
    }
  }
  boundary (h1);
  foreach_dimension()
    boundary (q1.x);
}

double dt = 0.;

void run (void)
{
  parameters();
  init_grid(N);

  events_init();
  init();
  foreach_dimension()
    boundary (q.x);
  foreach()
    h[] = max(h[], 0.);
  boundary (h);
  boundary (zb);

  timer_t start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    (* gradient) (zb, gzb); boundary (gzb.x); boundary (gzb.y);

    gradients();
    dt = dtnext (t, fluxes (DT));

    if (gradient == zero)
      /* 1st-order time-integration for 1st-order spatial
	 discretisation */
      update (q, q, h, h, dt);
    else {
      /* 2nd-order time-integration */
      /* predictor */
      update (q, q1, h, h1, dt/2.);
      
      vswap (q, q1);
      swap (h, h1);
      
      /* corrector */
      gradients();
      fluxes (dt);
      update (q1, q, h1, h, dt);
    }

    foreach(reduction(+:tnc)) tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
