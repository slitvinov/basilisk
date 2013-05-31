// 2D run-up of a solitary wave on a conical island, see e.g:
// [1] Hou, Liang, Simons and Hinkelman, 2013, section 3.4
//     http://dx.doi.org/10.1016/j.compfluid.2013.04.015
// [2] I.K. Nikolos, A.I. Delis, 2009, section 6.4.
//     http://dx.doi.org/10.1016/j.cma.2009.08.006

#include "saint-venant1.h"

#define MAXLEVEL 8
#define MINLEVEL 0

void parameters()
{
  L0 = 27.6;
  G = 9.81;
  N = 1 << MAXLEVEL;
}

// parameters for the solitary wave
// case B of [2]
double D = 0.32, H = 0.032, T = 2.45;
#define C sqrt(G*(D + H))

double sech2 (double x)
{
  double e = exp (-x);
  double s = 2.*e/(1. + e*e);
  return s*s;
}

double eta0 (double t)
{
  // eq (17) of [2]
  return H*sech2 (sqrt(3.*H/(4.*D*D*D))*C*(t - T));
}

double uleft (double t)
{
  double e = eta0 (t);
  return C*e/(D + e);
}

// wave paddle is on the left side
h[left] = D + eta0(t);
u.x[left] = uleft (t);
u.y[left] = 0.;

double island (double x, double y)
{
  // parameters for the conical island
  double r0 = 3.6, r1 = 1.1, H0 = 0.625;
  x -= 25.92/2.; y -= L0/2.;
  double r = sqrt(x*x + y*y);
  return r > r0 ? 0. : r < r1 ? H0 : (r - r0)/(r1 - r0)*H0;
}

#if QUADTREE
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = island (x, y);
}
#endif

// storage for maximum wave height
scalar hmax[];

void init()
{
#if QUADTREE
  zb.refine = refine_zb; // updates terrain
  // h refinement conserves elevation
  conserve_elevation();
#endif
  // initial conditions
  foreach() {
    zb[] = island (x, y);
    h[] = max(0., D - zb[]);
    hmax[] = 0.;
  }
}

// wave gauges

Gauge gauges[] = {
  {"WG3",   6.82, 13.05},
  {"WG6",   9.36, 13.80},
  {"WG9",  10.36, 13.80},
  {"WG16", 12.96, 11.22},
  {"WG22", 15.56, 13.80},
  {NULL}
};

event logfile (i++) {
  // stats on water depth
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
  assert (s.min >= 0.);

  // store hmax
  foreach()
    if (h[] > hmax[])
      hmax[] = h[];

  // output of point gauges
  output_gauges (gauges, {eta});
}

event outputfile (t = {9, 12, 13, 14, 20})
{
  static int nf = 0;
  printf ("file: conical-%d\n", nf);
  output_field ({h,zb,hmax}, N, stdout, true);
  printf ("file: level-%d\n", nf);
  scalar l[];
  foreach()
    l[] = level;
  output_field ({l}, N, stdout, false);
  nf++;
}

#if QUADTREE
event adapt (i++) {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  scalar w[];
  wavelet (eta, w);

  double cmax = 3e-4;
  int nf = refine_wavelet (w, cmax, MAXLEVEL, all);
  int nc = coarsen_wavelet (w, cmax/4., MINLEVEL, all);
  if (nf || nc)
    boundary (all);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}
#endif

int main () { run(); }
