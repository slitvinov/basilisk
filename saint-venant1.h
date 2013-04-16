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

  /* kinetic solver */
  if (hL > dry) {
    double ci = sqrt(G*hL/2.);
    double Mp = max(uL + ci*SQRT3, 0.);
    double Mm = max(uL - ci*SQRT3, 0.);
    double cig = ci/(6.*G*SQRT3);
    fh.x[i,0] = cig*3.*(Mp*Mp - Mm*Mm);
    fhu.x[i,0] = cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
    if (Mp > 0.) {
      double dt = CFL*DX/Mp;
      if (dt < dtmax)
	dtmax = dt;
    }
  }
  else
    fh.x[i,0] = fhu.x[i,0] = 0.;

  if (hR > dry) {
    double ci = sqrt(G*hR/2.);
    double Mp = min(uR + ci*SQRT3, 0.);
    double Mm = min(uR - ci*SQRT3, 0.);
    double cig = ci/(6.*G*SQRT3);
    fh.x[i,0] += cig*3.*(Mp*Mp - Mm*Mm);
    fhu.x[i,0] += cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
    if (Mm < 0.) {
      double dt = CFL*DX/-Mm;
      if (dt < dtmax)
	dtmax = dt;
    }
  }

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
    h1[] = h2[] - dt*(fh.x[] - fh.x[1,0])/DX;
    hu1[] = hu2[] - dt*(fhu.x[] - fhu.x[1,0] + Sb.x[])/DX;
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
