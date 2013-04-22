/* An implementation of:
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 * and
 *    [2] Kurganov, A., and Guergana, P. "A second-order
 *    well-balanced positivity preserving central-upwind scheme for
 *    the Saint-Venant system." Communications in Mathematical
 *    Sciences 5.1 (2007): 133-160.
 */

#include "utils.h"

scalar hu = new scalar, w = new scalar, B = new scalar;
scalar hu1 = new scalar, w1 = new scalar;
// gradients of hu, w
vector ghu = new vector, gw = new vector;
// fluxes
vector fhu = new vector, fw = new vector;

// Default parameters
// acceleration of gravity
double G = 1.;
// gradient
void (* gradient) (scalar *, vector *) = generalized_minmod;
// user-provided functions
void parameters (void);
void init       (void);

static void positivity()
{
  /* ensure positivity of reconstruction: see section 2.2 of [2] */
  foreach() {
    if (w[] + delta*gw.x[]/2. < (B[1,0] + B[])/2.)
      gw.x[] = (B[1,0] + B[] - 2.*w[])/delta; // (2.15) of [2]
    else if (w[] - delta*gw.x[]/2. < (B[-1,0] + B[])/2.)
      gw.x[] = (2.*w[] - B[-1,0] - B[])/delta; // (2.16) of [2]
  }
}

/* avoid division by zero */
static double hu_over_h (double hu, double h, double epsilon)
{
  double h4 = sq(h)*sq(h);
  return h*hu/sqrt((h4 + max(h4, epsilon))/2.); // (2.17) of [2]
}

static double flux (Point point, int i, double dtmax)
{
  VARIABLES;
  delta /= 2.;
  double hup = hu[i,0] - delta*ghu.x[i,0], hum = hu[i-1,0] + delta*ghu.x[i-1,0];
  double wp  = w[i,0] - delta*gw.x[i,0], wm = w[i-1,0] + delta*gw.x[i-1,0];
  double Bpm = (B[i,0] + B[i-1,0])/2.;
  double hp = wp - Bpm, hm = wm - Bpm;
  double epsilon = 1e-6;
  double up = hu_over_h (hup, hp, epsilon), um = hu_over_h (hum, hm, epsilon);
  hup = hp*up; hum = hm*um;
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double a = max(ap, -am);
  if (a > 0.) {
    fw.x[i,0] = (ap*hum - am*hup + ap*am*(wp - wm))/(ap - am); // (4.5) of [1]
    fhu.x[i,0] = (ap*(hum*um + G*sq(hm)/2.) - am*(hup*up + G*sq(hp)/2.) + 
		  ap*am*(hup - hum))/(ap - am);
    double dt = CFL*DX/a;
    return dt < dtmax ? dt : dtmax;
  }
  fw.x[i,0] = fhu.x[i,0] = 0.;
  return dtmax;
}

static void update (scalar hu2, scalar hu1, scalar w2, scalar w1, double dt)
{
  foreach() {
    double Bp = (B[1,0] + B[])/2., Bm = (B[-1,0] + B[])/2.;
    double S = G*(w[] - (Bp + Bm)/2.)*(Bp - Bm); // (2.10) of [2]
    hu1[] = hu2[] + dt*(fhu.x[] - fhu.x[1,0] - S)/DX;
  }
  foreach()
    w1[] = w2[] + dt*(fw.x[] - fw.x[1,0])/DX;
  boundary (w1);
  boundary (hu1);
}

double dt = 0.;

void run (void)
{
  parameters();
  init_grid(N);
  init();
  boundary (hu);
  boundary (B);
  foreach()
    w[] = max(w[], (B[-1,0] + 2.*B[] + B[1,0])/4.);
  boundary (w);

  timer_t start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    /* 2nd-order predictor-corrector */

    /* predictor */
    (* gradient) (scalars (hu, w), vectors (ghu, gw));
    boundary (ghu.x, gw.x);

    dt = DT;
    foreach()
      dt = flux (point, 0, dt);
    foreach_boundary(right)
      dt = flux (point, 1, dt);
    dt = dtnext (t, dt);

    update (hu, hu1, w, w1, dt/2.);
    swap (scalar, hu, hu1);
    swap (scalar, w, w1);

    /* corrector */
    (* gradient) (scalars (hu, w), vectors (ghu, gw));
    positivity();
    boundary (ghu.x, gw.x);

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
