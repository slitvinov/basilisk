/* An implementation of:
 *    Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */

#include "utils.h"
#include "events.h"

scalar hu = new scalar, w = new scalar, B = new scalar;
scalar hu1 = new scalar, w1 = new scalar;
// gradients of hu, w
vector ghu = new vector, gw = new vector;
// fluxes
vector fhu = new vector, fw = new vector;

#include "boundary.h"

// Default parameters
// acceleration of gravity
double G = 1.;
// gradient
void (* gradient) (const scalar, vector) = generalized_minmod;
// user-provided functions
void parameters (void);
void init       (void);

double flux (Point point, int i, double dtmax)
{
  VARIABLES;
  double hup = hu[i,0] - ghu.x[i,0]/2., hum = hu[i-1,0] + ghu.x[i-1,0]/2.;
  double wp  = w[i,0]  - gw.x[i,0]/2.,  wm  = w[i-1,0]  + gw.x[i-1,0]/2.;
  double Bpm  = (B[i,0] + B[i-1,0])/2.;
  double hp = wp - Bpm, hm = wm - Bpm;
  double up = hup/hp,  um = hum/hm;
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  fw.x[i,0] = (ap*hum - am*hup + (ap*am)*(wp - wm))/(ap - am); // (4.5) of [1]
  fhu.x[i,0] = (ap*(hum*hum/hm + sq(hm)/2.) - 
		am*(hup*hup/hp + sq(hp)/2.) + 
		(ap*am)*(hup - hum))/(ap - am);
  double a = max(ap, -am);
  double dt = CFL*L0*delta/a;
  if (dt < dtmax)
    fprintf (stderr, "%g %g %g %g %g\n", a, hp, hm, x, w[] - B[]);
  return dt < dtmax ? dt : dtmax;
}

#define swap(a,b) { scalar tmp = a; a = b; b = tmp; }

void update (scalar hu2, scalar hu1, scalar w2, scalar w1, double dt)
{
  foreach() {
    double Bp = (B[1,0] + B[])/2., Bm = (B[-1,0] + B[])/2.;
    double wp  = w[] + gw.x[]/2.,  wm  = w[] - gw.x[]/2.;
    double S = G*(Bp - Bm)*(wp - Bp + wm - Bm)/2.; // (4.4) of [1]
    hu1[] = hu2[] + dt*(fhu.x[] - fhu.x[1,0] - S)/(L0*delta);		  
  }
  foreach()
    w1[] = w2[] + dt*(fw.x[] - fw.x[1,0])/(L0*delta);
  boundary (w1);
  boundary (hu1);
}

double dt;

void run (void)
{
  parameters();
  init_grid(N);

  events_init();
  init();
  boundary (hu);
  boundary (w);
  boundary (B);

  timer_t start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    /* 2nd-order predictor-corrector */

    /* predictor */
    (* gradient) (hu, ghu); boundary (ghu.x);
    (* gradient) (w, gw);   boundary (gw.x);

    dt = DT;
    foreach()
      dt = flux (point, 0, dt);
    foreach_boundary(right)
      dt = flux (point, 1, dt);
    dt = dtnext (t, dt);

    update (hu, hu1, w, w1, dt/2.);
    swap (hu, hu1);
    swap (w, w1);

    /* corrector */
    (* gradient) (hu, ghu); boundary (ghu.x);
    (* gradient) (w, gw);   boundary (gw.x);

    foreach()
      flux (point, 0, dt);
    foreach_boundary(right)
      flux (point, 1, dt);

    update (hu1, hu, w1, w, dt);

    foreach() tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid ();
}
