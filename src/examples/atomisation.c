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

int maxlevel = 8;

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
  f.sigma = SIGMA;
  //  TOLERANCE = 1e-5;
  
  tt = timer_start();
  run();
}

event init (t = 0) {

  foreach()
    f[] = rhov[] = f0[] = 0.;
  foreach_face()
    alphav.x[] = muv.x[] = 0.;
  boundary ({f,f0,rhov,alphav,muv});
  
  refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel,
	  all);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(radius) - sq(y) - sq(z);
  fractions (phi, f0);
  
  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0}); // for boundary conditions on restricted f
 
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
}

/**
If gfsview is installed on the system, we can also visualise the
simulation as it proceeds. */

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen("gfsview2D ../atomisation.gfv", "w");
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
  scalar pid[], ff[];
  foreach() {
    pid[] = fmod(pid()*(npe() + 37), npe());
    ff[] = f[] < 1e-4 ? 0 : f[] > 1. - 1e-4 ? 1. : f[];
  }
  boundary ({pid,ff});
  output_gfs (file = name, t = t);
#endif
}

#ifdef DEBUGCOND
static void check_restriction (scalar a)
{
  double val = 123;
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  for (int l = 0; l <= maxlevel; l++) {
    if (l == 0)
      foreach_level_or_leaf (l)
	a[] = val;
    else
      foreach_level (l) {
	for (int i = -2; i <= 2; i++)
	  for (int j = -2; j <= 2; j++)
	    // fixme: boundary conditions should work too
	    assert ((aparent(i,j).pid < 0) || coarse(a,i,j) == val);
	a[] = val;
      }
    boundary_level ({a}, l);
  }
}

static void abortion (face vector uu)
{
  FILE * fp = lfopen ("uu", "w");
  foreach_face(x)
    for (int i = -2; i <= 2; i++)
      if (neighbor(0,i).pid >= 0)
	fprintf (fp, "%g %g %g %d\n", x, y + i*Delta, uu.x[0,i],
		 neighbor(0,i).pid);
  foreach_face(y)
    for (int i = -2; i <= 2; i++)
      if (neighbor(i).pid >= 0)
	fprintf (fp, "%g %g %g %d\n", x + i*Delta, y, uu.y[i], neighbor(i).pid);
  fclose (fp);
}
#endif // DEBUGCOND

event adapt (i++) {
  
#ifdef DEBUGCOND
  scalar ss[];
  face vector uu[];
  quadtree_trash ({ss,uu});
  foreach()
    ss[] = 1;
  foreach_face()
    uu.x[] = 1;
  boundary ({ss,uu});
#endif
  
  boundary ((scalar *){a});
  adapt_wavelet ({f,u}, (double[]){0.01,0.1,0.1,0.1}, maxlevel);
  restriction ({f0}); // for boundary conditions on restricted f
  event ("properties");
  
#ifdef DEBUGCOND
  foreach()
    foreach_neighbor()
      if (cell.pid >= 0) // fixme
	assert (ss[] == 1);
  check_restriction (ss);
  foreach_face()
    for (int i = -2; i <= 2; i++)
      if (neighbor(0,i).pid >= 0) {
	if (uu.x[0,i] != 1) {
	  abortion (uu);
	  fprintf (stderr, "uu.x = %g at %g %g\n", uu.x[0,i], x, y);
	}
	assert (uu.x[0,i] == 1);
      }
#endif
}
