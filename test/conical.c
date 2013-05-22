// 2D run-up of a solitary wave on a conical island
// see e.g. Hou, Liang, Simons and Hinkelman, 2013, section 3.4
// http://dx.doi.org/10.1016/j.compfluid.2013.04.015

#include "grid/cartesian.h"
#include "saint-venant1.h"

void parameters()
{
  L0 = 27.6;
  G = 9.81;
  N = 128;
}

double sech2 (double x)
{
  double e = exp (-x);
  double s = 2.*e/(1. + e*e);
  return s*s;
}

// parameters for the solitary wave
double D = 0.32, H = 0.064, T = 3.45;
#define C sqrt(G*(D + H))

double eta (double t)
{
  return H*sech2 (sqrt(3.*H/(4.*D))*C*(t - T));
}

double uleft (double t)
{
  double e = eta (t);
  return C*e/(D + e);
}

h[left] = D + eta(t);
u.x[left] = uleft (t);
u.y[left] = 0.;

scalar hmax[];

// wave gauges
FILE * WG3, * WG6, * WG9, * WG16, * WG22;

void init()
{
  // parameters for the conical island
  double r0 = 3.6, r1 = 1.1, H0 = 0.625;
  foreach() {
    x -= 25.92/2.; y -= L0/2.;
    double r = sqrt(x*x + y*y);
    zb[] = r > r0 ? 0. : r < r1 ? H0 : (r - r0)/(r1 - r0)*H0;
    h[] = max(0., D - zb[]);
    hmax[] = 0.;
  }
  WG3 = fopen ("WG3", "w");
  WG6 = fopen ("WG6", "w");
  WG9 = fopen ("WG9", "w");
  WG16 = fopen ("WG16", "w");
  WG22 = fopen ("WG22", "w");
}

int event (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);

  foreach()
    if (h[] > hmax[])
      hmax[] = h[];

  double x = 6.82, y = 13.05;
  fprintf (WG3, "%g %g\n", t,
	   interpolate (zb, x, y) + interpolate (h, x, y));
  x = 9.36, y = 13.80;
  fprintf (WG6, "%g %g\n", t,
	   interpolate (zb, x, y) + interpolate (h, x, y));
  x = 10.36, y = 13.80;
  fprintf (WG9, "%g %g\n", t,
	   interpolate (zb, x, y) + interpolate (h, x, y));
  x = 12.96, y = 11.22;
  fprintf (WG16, "%g %g\n", t,
	   interpolate (zb, x, y) + interpolate (h, x, y));
  x = 15.56, y = 13.80;
  fprintf (WG22, "%g %g\n", t,
	   interpolate (zb, x, y) + interpolate (h, x, y));
}

int event (t = {9, 11, 13, 20})
{
  static int nf = 0;
  printf ("file: conical-%d\n", nf++);
  output_field ({h,zb,hmax}, N, stdout, true);
}

int main () { run(); }
