// #include "grid/multigrid.h"
//#include "navier-stokes/mac.h"
#include "navier-stokes/centered.h"
// #include "bcg.h"
#include "vof.h"

#define LEVEL 7

scalar f[];
scalar * interfaces = {f};

void parameters() {
  L0 = 4;
  // coordinates of lower-left corner
  X0 = Y0 = -2;
  // number of grid points
  N = 1 << LEVEL;
  // viscosity
  const face vector muc[] = {0.00313,0.00313};
  mu = muc;
  // maximum timestep
  DT = 5e-3;
}

event init (t = 0) {
  alpha = new face vector;
  scalar phi[];
  foreach_vertex()
    phi[] = 0.05*cos (2.*pi*x) + y;
  fractions (phi, f);
}

event advance (i++) {
  trash ({alpha});
  foreach_face() {
    double fm = (f[] + f[-1,0])/2.;
    alpha.x[] = 1./(fm*1.225 + (1. - fm)*0.1694);
  }
  boundary_normal ({alpha});
}

event projection (i++) {
  foreach()
    u.y[] -= 9.81*dt;
  boundary ((scalar *){u});
}

event logfile (i++) {
  stats s = statsf (f);
  printf ("%g %d %g %g %g %g %d %d %d\n", 
	  t, i, dt, s.sum - 8., s.min, s.max - 1., mgp.i, mgpf.i, mgu.i);
  assert (s.min >= -1e-10 && s.max <= 1. + 1e-10);
  assert (fabs (s.sum - 8.) < 1e-4);
}

event interface (t = {0,0.7,0.8,0.9,1.}) {
  char s[80];
  sprintf (s, "interface-%g", t);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);

  if (t == 1.)
    output_facets (f, stderr);

  sprintf (s, "field-%g", t);
  fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf}, fp);
  fclose (fp);
}

event movie (i += 3) {
  static FILE * fp = popen ("ppm2mpeg > f.mpg", "w");
  output_ppm (f, fp, spread = 2);

  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL);
}

void output_all (scalar * l, FILE * fp)
{
  foreach() {
    fprintf (fp, "%g %g", x, y);
    for (scalar s in l)
      fprintf (fp, " %g", s[]);
    fputc ('\n', fp);
  }
}

#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-3,0.1,0.1}, LEVEL, list = {p,u,pf,uf,f});
}
#endif

int main() { run (); }
