// #include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

#define radius 1./12.
#define length 0.025
#define Re 5800

#define rho1 1.

#define rho2 1./27.84
#define SIGMA 3e-5

/**
The density and viscosity are defined using the arithmetic average. */

#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f)  (2.*radius/Re*rho(f))

int maxlevel = 12;

scalar f[], * interfaces = {f};
face vector alphav[];
scalar rhov[];
face vector muv[];

scalar f0[];

u.n[left] = dirichlet(f0[]*(1. + 0.05*sin (10.*2.*pi*t)));
u.t[left] = dirichlet(0);
p[left]   = neumann(0);
f[left]   = f0[];

u.n[right] = neumann(0);
p[right] = dirichlet(0);

#if 1
p[top] = neumann(0);
uf.n[top]    = 0;
p[bottom] = neumann(0);
uf.n[bottom] = 0;

#if dimension == 3
p[front] = neumann(0);
uf.n[front]  = 0;
p[back] = neumann(0);
uf.n[back]   = 0;
#endif

#endif

timer tt;

int main (int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  
  init_grid (64);

  origin (0, -1.5, -1.5);
  size (3.);

  /**
  The density and viscosity are defined by the variable fields we
  allocated above. We also set the surface tension for interface
  *f*. */
  
  alpha = alphav;
  rho = rhov;
  mu = muv;
  //  f.sigma = SIGMA;
  TOLERANCE = 1e-5;
  
  tt = timer_start();
  run();
}

event init (t = 0) {

  refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel,
	  all);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(radius) - sq(y) - sq(z);
  fractions (phi, f0);
  f0.refine = f0.prolongation = fraction_refine;  

  foreach() {
    f[] = f0[]*(x < length);
    u.x[] = f[];
  }
  boundary ({f,u.x});
}

event properties (i++) {
#if 1
  f.prolongation = refine_bilinear;
  boundary ({f});
#endif
  
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);

#if 1
  f.prolongation = fraction_refine;
  boundary ({f});
#endif
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume. */

event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  int nc = 0;
  static long tnc = 0;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb) reduction(+:nc)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
    nc++;
  }
  double elapsed = timer_elapsed (tt);
  tnc += nc;
  fprintf (ferr, "%g %g %g %g %g %g %d %d %d %d %g %g\n", 
	   t, sb, -1., xb/sb, vb/sb, dt, mgp.i, mgpf.i, mgu.i,
	   nc, elapsed, tnc/elapsed);

#if DEBUG_MPI
  if (i >= 35)
    debug_iteration = i;
#endif
}

/**
If gfsview is installed on the system, we can also visualise the
simulation as it proceeds. */

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview3D ../atomisation.gfv", "w");
  output_gfs (fp, t = t, list = {u,p,f});
}
#endif

#if 0 //!_MPI
event movie (t += 4e-3)
{
  static FILE * fp = popen ("gfsview-batch3D ../atomisation.gfv "
			    "| ppm2mpeg -s 800x600 > atomisation.mpg", "w");
  output_gfs (fp, t = t, list = {u,p,f});
  fprintf (fp, "Save stdout { format = PPM width = 1600 height = 1200 }\n");
}
#endif

event snapshot (t = 0.1; t += 0.1; t <= 1.8) { // t <= 1.8
#if 1
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  scalar pid[];
  foreach()
    pid[] = pid();
  output_gfs (file = name, t = t);
#endif
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,0.1,0.1,0.1}, maxlevel);
  event ("properties");
}
