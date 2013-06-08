#include "saint-venant.h"
#include "terrain.h"
#include "okada.h"

#define MAXLEVEL 10
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HMAXE    5e-2 // error on maximum free surface elevation (5 cm)

// metres to degrees
double mtd = 360./40075e3;

void parameters()
{
#if QUADTREE
  // 32^2 grid points to start with
  N = 1 << MINLEVEL;
#else // Cartesian
  // 1024^2 grid points
  N = 1 << MAXLEVEL;
#endif
  // the domain is 54 degrees squared
  L0 = 54.;
  // centered on 94,8 longitude,latitude
  X0 = 94 - L0/2.;
  Y0 = 8. - L0/2.;
  /* rescale G so that time is in minutes, horizontal length scales in
     degrees and vertical length scales in metres */
  G = 9.81*sq(mtd)*sq(60.);
}

// "radiation" boundary conditions on left,right,bottom
// sealevel at zero
u.x[left]   = - radiation(0);
u.x[right]  = + radiation(0);
u.y[bottom] = - radiation(0);

// extra storage for hmax (maximum wave elevation)
scalar hmax[];

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  // adapt on both eta and hmax with errors ETAE and HMAXE respectively
  astats s = adapt_wavelet ({eta, hmax}, (double[]){ETAE, HMAXE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else
  return 0;
#endif
}

void init()
{
  // use etopo2 kdt terrain database for topography zb
  terrain (zb, "/home/popinet/terrain/etopo2", NULL);
  // ensure that water level is conserved during refinement/coarsening
  // the default is to conserve volume
  conserve_elevation();
  // initial still water surface is at z = 0 i.e. the water depth h[] is...
  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});
  // initial deformation
  fault (x = 94.57, y = 3.83,
	 depth = 11.4857e3,
	 strike = 323, dip = 12, rake = 90,
	 length = 220e3, width = 130e3,
	 U = 18,
	 iterate = adapt);
}

// second fault segment is triggered at t = 272 seconds
event fault2 (t = 272./60.) {
  fault (x = 93.90, y = 5.22,
	 depth = 11.4857e3,
	 strike = 348, dip = 12, rake = 90,
	 length = 150e3, width = 130e3,
	 U = 23,
	 iterate = adapt);
}

// third fault segment is triggered at t = 588 seconds
event fault3 (t = 588./60.)
{
  fault (x = 93.21, y = 7.41,
	 depth = 12.525e3,
	 strike = 338, dip = 12, rake = 90,
	 length = 390e3, width = 120e3,
	 U = 12,
	 iterate = adapt);
}

// fourth fault segment is triggered at t = 913 seconds
event fault4 (t = 913./60.)
{
  fault (x = 92.60, y = 9.70,
	 depth = 15.12419e3,
	 strike = 356, dip = 12, rake = 90,
	 length = 150e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
}

// fifth fault segment is triggered at t = 1273 seconds
event fault5 (t = 1273./60.)
{
  fault (x = 92.87, y = 11.70,
	 depth = 15.12419e3,
	 strike = 10, dip = 12, rake = 90,
	 length = 350e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
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
    // hmax
    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
  }
  boundary ({hmax, u});
}

// snapshots every hour
event snapshots (t += 60; t <= 600) {
  printf ("file: t-%g\n", t);
  output_field ({h, zb, hmax}, stdout, n = 1 << MAXLEVEL, linear = true);
}

// movies every minute
event movies (t++) {
  static FILE * fp = NULL;
  if (!fp) fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, mask = m, min = -2, max = 2, n = 512, linear = true);

  static FILE * fp2 = NULL;
  if (!fp2) fp2 = popen ("ppm2mpeg > eta-zoom.mpg", "w");
  output_ppm (etam, fp2, mask = m, min = -2, max = 2, n = 512, linear = true,
	      box = {{89,8},{98,16}});

  static FILE * fp1 = NULL;
  if (!fp1) fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = MAXLEVEL, n = 512);

#if _OPENMP
  static FILE * fp3 = NULL;
  if (!fp3) fp3 = popen ("ppm2mpeg > pid.mpg", "w");
  foreach()
    etam[] = pid();
  double tmax = omp_get_max_threads() - 1;
  output_ppm (etam, fp3, max = tmax, n = 512);
#endif // _OPENMP
}

// tide gauges

Gauge gauges[] = {
  // file   lon      lat         description
  {"coco", 96.88,  -12.13, "Cocos Islands, Australia"},
  {"colo", 79.83,    6.93, "Colombo, Sri Lanka"},
  {"hani", 73.18,    6.77, "Hanimaadhoo, Maldives"},
  {"male", 73.52,    4.18, "Male, Maldives"},
  {"gana", 73.17,   -0.68, "Gan, Maldives"},
  {"dieg", 72.38,    -7.3, "Diego Garcia, UK"},
  {"rodr", 63.42,  -19.67, "Rodriguez I., Mauritius"},
  {"loui", 57.5,   -20.15, "Port Louis, Mauritius"},
  {"lare", 55.3,   -20.92, "La Reunion, France"},
  {"hill", 115.73, -31.82, "Hillarys, Australia"},
  {"sala", 54,         17, "Salalah, Oman"},
  {"laru", 55.53,   -4.68, "Pointe La Rue, Seychelles"},
  {"lamu", 40.9,    -2.27, "Lamu, Kenya"},
  {"zanz", 39.18,   -6.15, "Zanzibar, Tanzania"},
  {"chen", 80.3,     13.1, "Chennai, India"},
  {"para", 86.7,    20.26, "Paradip, India"},
  {"visa", 83.28,   17.68, "Visakhapatnam, India"},
  {"koch", 76.26,    9.96, "Kochi, India"},
  {"morm", 73.8,    15.42, "Mormugao, India"},
  {"okha", 69.08,   22.47, "Okha, India"},
  {"tuti", 78.15,     8.8, "Tuticorin, India"},
  {"taru", 99.65,   6.702, "Tarutao, Thailand"},
  {"tapa", 98.425,  7.765, "Tapaonoi, Thailand"},
  {NULL}
};

event gauges1 (i++) output_gauges (gauges, {eta});

event do_adapt (i++) adapt();

int main() { run(); }
