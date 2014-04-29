#include "saint-venant.h"
#include "terrain.h"

#define MAXLEVEL 9
#define MINLEVEL 5

// metres to degrees
double mtd = 360./40075e3;

int main run()
{
  // 512^2 grid points
  init_grid (1 << MAXLEVEL);
  // the domain is 5 degrees squared
  size (5.);
  // centered on 174,-40.8 longitude,latitude
  origin (174 - L0/2., -40.8 - L0/2.);
  /* rescale G so that time is in minutes, horizontal length scales in
     degrees and vertical length scales in metres */
  G = 9.81*sq(mtd)*sq(60.);
}

// M2 tidal frequency. The period is 12h25 (745 minutes).
#define M2F (2.*M_PI/745.)

// "radiation" boundary conditions on left,right,top,bottom
u.x[left]   = - radiation (  sin(M2F*t));
u.x[top]    =   radiation (  sin(M2F*t));
u.x[right]  = + radiation (- sin(M2F*t));
u.x[bottom] = - radiation (- sin(M2F*t));

void init()
{
  // use several databases for topography zb
  terrain (zb, 
	   "/home/popinet/terrain/niwa", 
	   "/home/popinet/terrain/geo",
	   "/home/popinet/terrain/etopo2", 
	   NULL);
  // ensure that water level is conserved during refinement/coarsening
  // the default is to conserve volume
  conserve_elevation();
  // sealevel at z = 0
  foreach()
    h[] = max(0., - zb[]);
}

// every timestep
event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);

  foreach() {
    // quadratic bottom friction, coefficient 1e-4 (dimensionless)
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/(h[]*mtd);
    foreach_dimension()
      u.x[] /= a;
  }
  boundary ((scalar *){u});
}

// snapshots every hour
event snapshots (t += 60; t <= 6000) {
  printf ("file: t-%g\n", t);
  output_field ({h, zb, eta, u}, stdout, n = N, linear = true);
}

// movies every 10 minutes
event movies (t += 10) {
  static FILE * fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, mask = m, min = -1, max = 1, n = 512, linear = true);

  static FILE * fp2 = popen ("ppm2mpeg > vort.mpg", "w");
  scalar vort = etam;
  foreach()
    vort[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*delta);
  boundary ({vort});
  output_ppm (vort, fp2, mask = m, // min = -1e-2, max = 1e-2, 
	      n = 512, linear = true);

  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = MAXLEVEL, n = 512);
}

// tide gauges

Gauge gauges[] = {
  {"well", 96.88,  -12.13, "Wellington, New Zealand"},
  {NULL}
};

event gauges1 (i++) output_gauges (gauges, {eta});

event adapt (i++) {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  scalar w[];
  wavelet (eta, w);

  double cmax = 1e-2;
  int nf = refine_wavelet (w, cmax, MAXLEVEL, all);
  int nc = coarsen_wavelet (w, cmax/4., MINLEVEL, all);
  if (nf || nc)
    boundary (all);

  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", nf, nc);
}
