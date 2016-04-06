#define JACOBI 1 // fixme: does not converge without this
#include "navier-stokes/centered.h"
#include "vof.h"

#define LEVEL 8
#define BGHOSTS 2 // fixme: need this to avoid FPEs when using -catch

scalar f[], * interfaces = {f};
face vector alphav[];
scalar rhov[];

uf.n[left]   = 0;
p[left] = neumann(0.);
uf.n[right]  = 0;
p[right] = neumann(0.);

int main() {
  size (4);
  origin (-2, -2);
  init_grid (1 << LEVEL);
  // viscosity
  const face vector muc[] = {0.00313,0.00313};
  mu = muc;
  DT = 5e-3;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
  mask (x < -0.5 ? left : x > 0.5 ? right : none);
  
  const face vector g[] = {0,-9.81};
  a = g;
  alpha = alphav;
  rho = rhov;

  vertex scalar phi[];
  foreach_vertex()
    phi[] = 0.05*cos (2.*pi*x) + y;
  fractions (phi, f);
}

#define rho(f) ((f)*1.225 + (1. - (f))*0.1694)

event properties (i++) {
  trash ({alphav,rhov});
  foreach_face() {
    double ff = (f[] + f[-1,0])/2.;
    alphav.x[] = fm.x[]/rho(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);
}

event logfile (i++) {
  stats s = statsf (f);
  printf ("%g %d %g %g %g %g %d %d %d\n", 
	  t, i, dt, s.sum - 8., s.min, s.max - 1., mgp.i, mgpf.i, mgu.i);
  assert (s.min >= -1e-10 && s.max <= 1. + 1e-10);
  assert (fabs (s.sum - 2.) < 5e-5);
}

// event interface (t = {0,0.7,0.8,0.9,1.}) {
event interface (t = 0.7) {
  output_facets (f, stderr);
  FILE * fp = fopen ("levels", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_field ({l}, fp);
  fclose (fp);
}

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview2D -s rt.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){5e-3}, LEVEL);
}
#endif
