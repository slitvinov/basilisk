#include "utils.h"
#include "events.h"

scalar hu = new scalar, h = new scalar, zb = new scalar;
scalar hu1 = new scalar, h1 = new scalar;
// topographic source term
vector Sb = new vector;
// gradients of hu, h, zb
vector ghu = new vector, gh = new vector, gzb = new vector;
// fluxes
vector fhu = new vector, fh = new vector;

#include "boundary.h"

// Default parameters
// acceleration of gravity
double G = 1.;
// gradient
void (* gradient) (const scalar, vector) = generalized_minmod;
// dry
double dry = 1e-6;
// user-provided functions
void parameters (void);
void init       (void);

#define SQRT3 1.73205080756888

/* kinetic(), kurganov() and hllc() are all Riemann solvers */

void kinetic (double hL, double hR, double uL, double uR, double delta,
	      double * fh, double * fhu, double * dtmax)
{
  double ci = sqrt(G*hL/2.);
  double Mp = max(uL + ci*SQRT3, 0.);
  double Mm = max(uL - ci*SQRT3, 0.);
  double cig = ci/(6.*G*SQRT3);
  *fh = cig*3.*(Mp*Mp - Mm*Mm);
  *fhu = cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
  if (Mp > 0.) {
    double dt = CFL*delta/Mp;
    if (dt < *dtmax)
      *dtmax = dt;
  }

  ci = sqrt(G*hR/2.);
  Mp = min(uR + ci*SQRT3, 0.);
  Mm = min(uR - ci*SQRT3, 0.);
  cig = ci/(6.*G*SQRT3);
  *fh += cig*3.*(Mp*Mp - Mm*Mm);
  *fhu += cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
  if (Mm < 0.) {
    double dt = CFL*delta/-Mm;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

void kurganov (double hm, double hp, double um, double up, double delta,
	       double * fh, double * fhu, double * dtmax)
{
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double hum = hm*um, hup = hp*up;
  double a = max(ap, -am);
  if (a > 0.) {
    *fh = (ap*hum - am*hup + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fhu = (ap*(hum*um + G*sq(hm)/2.) - am*(hup*up + G*sq(hp)/2.) + 
	    ap*am*(hup - hum))/(ap - am);
    double dt = CFL*delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fhu = 0.;
}

void hllc (double hL, double hR, double uL, double uR, double delta,
	   double * fh, double * fhu, double * dtmax)
{
  double cL = sqrt (G*hL), cR = sqrt (G*hR);
  double ustar = (uL + uR)/2. + cL - cR;
  double cstar = (cL + cR)/2. + (uL - uR)/4.;
  double SL = hL == 0. ? uR - 2.*cR : min (uL - cL, ustar - cstar);
  double SR = hR == 0. ? uL + 2.*cL : max (uR + cR, ustar + cstar);

  if (0. <= SL) {
    *fh = uL*hL;
    *fhu = hL*(uL*uL + G*hL/2.);
  }
  else if (0. >= SR) {
    *fh = uR*hR;
    *fhu = hR*(uR*uR + G*hR/2.);
  }
  else {
    double fhL = uL*hL;
    double fhuL = hL*(uL*uL + G*hL/2.);
    double fhR = uR*hR;
    double fhuR = hR*(uR*uR + G*hR/2.);
    *fh = (SR*fhL - SL*fhR + SL*SR*(hR - hL))/(SR - SL);
    *fhu = (SR*fhuL - SL*fhuR + SL*SR*(hR*uR - hL*uL))/(SR - SL);
#if 0
    double SM = ((SL*hR*(uR - SR) - SR*hL*(uL - SL))/
		  (hR*(uR - SR) - hL*(uL - SL)));
    if (SL <= 0. && 0. <= SM)
      f[V] = uL[V]*f[H];
    else if (SM <= 0. && 0. <= SR)
      f[V] = uR[V]*f[H];
    else
      assert (false);
#endif
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > 0.) {
    double dt = CFL*delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

static double flux (Point point, int i, double dtmax)
{
  VARIABLES;
  double eta = h[i,0], etan = h[i-1,0];
  if (eta <= dry && etan <= dry) {
    fh.x[i,0] = fhu.x[i,0] = 0.;
    return dtmax;
  }
  
  double zbL = zb[i,0] - gzb.x[i,0]/2.;
  double zbR = zb[i-1,0] + gzb.x[i-1,0]/2.;
  double zbLR = max(zbL, zbR);

  double etaL = eta <= dry ? 0. : eta - gh.x[i,0]/2.;
  double uL = etaL <= dry ? 0. : (hu[i,0] - ghu.x[i,0]/2.)/etaL;
  double hL = max(0., etaL + zbL - zbLR);
    
  double etaR = etan <= dry ? 0. : etan + gh.x[i-1,0]/2.;
  double uR = etaR <= dry ? 0. : (hu[i-1,0] + ghu.x[i-1,0]/2.)/etaR;
  double hR = max(0., etaR + zbR - zbLR);

  kurganov (hR, hL, uR, uL, DX, &fh.x[i,0], &fhu.x[i,0], &dtmax);

  /* topographic source term */
  if (eta <= dry) eta = 0.;
  if (etan <= dry) etan = 0.;
  Sb.x[i,0] -= G/2.*(sq(hL) - sq(etaL) + (etaL + eta)*gzb.x[i,0]/2.);
  Sb.x[i-1,0] += G/2.*(sq(hR) - sq(etaR) - (etaR + etan)*gzb.x[i-1,0]/2.);

  return dtmax;
}

#define swap(a,b) { scalar tmp = a; a = b; b = tmp; }

static void update (scalar hu2, scalar hu1, scalar h2, scalar h1, double dt)
{
  foreach() {
    h1[] = h2[] + dt*(fh.x[] - fh.x[1,0])/DX;
    hu1[] = hu2[] + dt*(fhu.x[] - fhu.x[1,0] + Sb.x[])/DX;
    Sb.x[] = 0.;
  }
  boundary (h1);
  boundary (hu1);
}

double dt = 0.;

void run (void)
{
  parameters();
  init_grid(N);

  events_init();
  init();
  boundary (hu);
  boundary (h);
  boundary (zb);

  timer_t start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    (* gradient) (zb, gzb); boundary (gzb.x);

    /* 2nd-order predictor-corrector */
    
    /* predictor */
    (* gradient) (hu, ghu); boundary (ghu.x);
    (* gradient) (h, gh); boundary (gh.x);

    dt = DT;
    foreach()
      dt = flux (point, 0, dt);
    foreach_boundary(right)
      dt = flux (point, 1, dt);
    dt = dtnext (t, dt);

    update (hu, hu1, h, h1, dt/2.);
    swap (hu, hu1);
    swap (h, h1);

    /* corrector */
    (* gradient) (hu, ghu); boundary (ghu.x);
    (* gradient) (h, gh); boundary (gh.x);

    foreach()
      flux (point, 0, dt);
    foreach_boundary(right)
      flux (point, 1, dt);

    update (hu1, hu, h1, h, dt);

    foreach() tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
