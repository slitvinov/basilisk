#include "grid/multigrid.h"
#include "saint-venant1.h"
#include "terrain.h"
#include "okada.h"

// metres to degrees
double mtd = 360./40075e3;

void parameters()
{
  // 1024^2 grid points
  N = 1024;
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
u.x[left]   = - 2.*(sqrt (G*h[]) - sqrt(G*max(-zb[], 0.)));
u.x[right]  = + 2.*(sqrt (G*h[]) - sqrt(G*max(-zb[], 0.)));
u.y[bottom] = - 2.*(sqrt (G*h[]) - sqrt(G*max(-zb[], 0.)));

// extra storage for hmax (maximum wave elevation)
scalar hmax[];

void init()
{
  // use etopo2 kdt terrain database for topography zb
  terrain (zb, "/home/popinet/terrain/etopo2");
  //  zb.refine = elevation;
  scalar d[];
  foreach()
    d[] = hmax[] = 0.;
  // initial deformation
  okada (d, 
	 .x = 94.57, .y = 3.83,
	 .depth = 11.4857e3,
	 .strike = 323, .dip = 12, .rake = 90,
	 .length = 220e3, .width = 130e3,
	 .U = 18);
  // sealevel at z = 0 + initial deformation
  foreach() {
    h[] = max(0., - zb[]);
    if (h[] > dry)
      h[] = max (0., h[] + d[]);
  }    
}

// second fault segment is triggered at t = 272 seconds
int event (t = 272./60.)
{
  scalar d[];
  foreach()
    d[] = 0.;
  okada (d,
	 .x = 93.90, .y = 5.22,
	 .depth = 11.4857e3,
	 .strike = 348, .dip = 12, .rake = 90,
	 .length = 150e3, .width = 130e3,
	 .U = 23);
  //deformation is added to h[] (water depth) only in wet areas
  foreach()
    if (h[] > dry)
      h[] = max (0., h[] + d[]);
}

// third fault segment is triggered at t = 588 seconds
int event (t = 588./60.)
{
  scalar d[];
  foreach()
    d[] = 0.;
  okada (d,
	 .x = 93.21, .y = 7.41,
	 .depth = 12.525e3,
	 .strike = 338, .dip = 12, .rake = 90,
	 .length = 390e3, .width = 120e3,
	 .U = 12);
  foreach()
    if (h[] > dry)
      h[] = max (0., h[] + d[]);
}

// fourth fault segment is triggered at t = 913 seconds
int event (t = 913./60.)
{
  scalar d[];
  foreach()
    d[] = 0.;
  okada (d,
	 .x = 92.60, .y = 9.70,
	 .depth = 15.12419e3,
	 .strike = 356, .dip = 12, .rake = 90,
	 .length = 150e3, .width = 95e3,
	 .U = 12);
  foreach()
    if (h[] > dry)
      h[] = max (0., h[] + d[]);
}

// fifth fault segment is triggered at t = 1273 seconds
int event (t = 1273./60.)
{
  scalar d[];
  foreach()
    d[] = 0.;
  okada (d,
	 .x = 92.87, .y = 11.70,
	 .depth = 15.12419e3,
	 .strike = 10, .dip = 12, .rake = 90,
	 .length = 350e3, .width = 95e3,
	 .U = 12);
  foreach()
    if (h[] > dry)
      h[] = max (0., h[] + d[]);
}

// at every timestep
int event (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %.8f %g %g %g\n", t, i, s.min, s.max, s.sum, 
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
  boundary ((scalar *){u});
}

// snapshots every hour
int event (t += 60; t <= 600) {
  static int nf = 0;
  printf ("file: file-%d\n", nf++);
  output_field ({h, zb, hmax}, N, stdout, true);
}

// movie every minute
int event (t++) {
  static FILE * ppm = NULL;
  if (!ppm) ppm = popen ("ppm2mpeg > eta.mpg", "w");
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : nodata;
  boundary ({eta});
  output_ppm (eta, -2, 2, 512, ppm, false);
}

int main() { run(); }
