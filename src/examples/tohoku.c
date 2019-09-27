/**
# The 2004 Indian Ocean tsunami

The 2004 Indian Ocean tsunami was caused by a large-scale fault
rupture (> 1000 km) at the Indian–Australian and Eurasian–Andaman
plate boundaries. This example uses the fault model of [Grilli et
al, 2007](/src/references.bib#grilli2007) as initial conditions for a 
Saint-Venant solution of the subsequent tsunami. A similar setup is
discussed in [Popinet, 2011](/src/references.bib#popinet2011).

## Solver setup

The following headers specify that we use the [Saint-Venant
solver](/src/saint-venant.h) together with [(dynamic) terrain
reconstruction](/src/terrain.h) and the [Okada fault
model](/src/okada.h). */

#include "spherical.h"
#if ML
# include "layered/hydro.h"
# include "layered/nh-box.h"
scalar h;
vector u;
# define MGD mgp.i
#else
#if 1
# include "green-naghdi.h"
# define MGD mgD.i
#else
# include "saint-venant.h"
#endif
#endif
#include "terrain.h"
#include "okada.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 13
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HMAXE    5e-2 // error on maximum free surface elevation (5 cm)

int main()
{
  /**
  Here we setup the domain geometry. For the moment Basilisk only
  supports square domains. For this example we have to use degrees as
  horizontal units because that is what the topographic database uses
  (eventually coordinate mappings will give more flexibility). We set
  the size of the box *L0* and the coordinates of the lower-left corner
  *(X0,Y0)*. */

  Radius = 6371220.;
  // the domain is 73 degrees squared
  size (73.);
  // centered on 142,38 longitude,latitude
  origin (142. - L0/2., 38. - L0/2.);

  /**
  *G* is the acceleration of gravity required by the Saint-Venant
  solver. This is the only dimensional parameter. We rescale it so that
  time is in minutes, horizontal distances in degrees and vertical
  distances in metres. This is a trick to circumvent the current lack of
  coordinate mapping. */

  // acceleration of gravity in degrees^2/min^2/m
  G = 9.81*sq(60.);

  /**
  When using a quadtree (i.e. adaptive) discretisation, we want to start
  with the coarsest grid, otherwise we directly refine to the maximum
  level. Note that *1 << n* is C for $2^n$. */

#if QUADTREE
  // 32^2 grid points to start with
  init_grid (1 << MINLEVEL);
#else // Cartesian
  // 1024^2 grid points
  init_grid (1 << MAXLEVEL);
#endif

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

#if ML
  filter = 3.; // 3 minutes filter
#if NH
  breaking = 0.07;
  NITERMIN = 1;
  TOLERANCE = HUGE;
#endif
#endif

  run();
}

/**
We declare and allocate another scalar field which will be used to
store the maximum wave elevation reached over time. */

scalar hmax[];

/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two *#if...#else* branches selecting
whether the simulation is being run on an (adaptive) quadtree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height *hmax* -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
*η*) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int my_adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  /**
  We can now use wavelet adaptation on the list of scalars *{η,hmax}*
  with thresholds *{ETAE,HMAXE}*. The compiler is not clever enough yet
  and needs to be told explicitly that this is a list of *double*s,
  hence the *(double[])*
  [type casting](http://en.wikipedia.org/wiki/Type_conversion). 
  
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({eta, hmax}, (double[]){ETAE,HMAXE},
			    MAXLEVEL, MINLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

/**
## Initial conditions

We first specify the terrain database to use to reconstruct the
topography $z_b$. This KDT database needs to be built beforehand. See the
[*xyz2kdt* manual](http://gfs.sourceforge.net/wiki/index.php/Xyz2kdt)
for explanations on how to do this.

The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */

event init (i = 0)
{
#if ML
  h = hl[0];
  u = ql[0];
#endif
  
  terrain (zb, "~/terrain/etopo2", "~/terrain/srtm_japan", NULL);
  conserve_elevation();

  /**
  The initial still water surface is at $z=0$ so that the water depth
  $h$ is... */

  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});

  /**
  The initial deformation is given by an Okada fault model with the
  following parameters. The *iterate = adapt* option will iterate this
  initialisation until our *adapt()* function above returns zero
  i.e. until the deformations are resolved properly. */
  
  #include "faults1.h"

  /**
  ## Boundary conditions

  We set the normal velocity component on the left, right and bottom
  boundaries to a "radiation condition" with a reference sealevel of
  zero. The top boundary is always "dry" in this example so can be left
  alone. Note that the sign is important and needs to reflect the
  orientation of the boundary. */

  u.n[left]   = - radiation(0);
  u.n[right]  = + radiation(0);
  u.t[bottom] = - radiation(0);
  u.t[top]    = + radiation(0);
}

/**
## Outputs

### At each timestep

We output simple summary statistics for *h* and *u.x* on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (ferr,
	     "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i speed tn\n");
  fprintf (ferr, "%g %d %g %g %g %g %g %g %d %g %ld\n",
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD,
	   perf.speed, grid->tn);
}

event viscous_term (i++)
{
  if (filter) {
#if 0
    assert (nl == 1); // fixme
    scalar h = hl[0];
    vector q = ql[0];
    foreach() {
      for (scalar s in {q})
	foreach_dimension()
	  if (h[-1] > dry && h[] > dry && h[1] > dry) {
	    double Dp = s[1]/h[1] - s[]/h[], Dm = s[]/h[] - s[-1]/h[-1];
	    if (Dp*Dm < 0. && ((s[2]/(h[2] + dry) + s[]/h[] - 2.*s[1]/h[1])*
			       (s[-1]/h[-1] + s[1]/h[1] - 2.*s[]/h[]) < 0. ||
			       (s[-1]/h[-1] + s[1]/h[1] - 2.*s[]/h[])*
			       (s[-2]/(h[-2] + dry) + s[]/h[] - 2.*s[-1]/h[-1]) < 0.)) {
	      double dp, dm;
	      if (fabs(Dp) > fabs(Dm)) {
		dp = fabs(Dp);
		dm = fabs(Dm);
	      }
	      else {
		dp = fabs(Dm);
		dm = fabs(Dp);
	      }
	      double d = min(dm, dp/2.);
	      double a = Dp > 0. ? 1. : -1.;
	      s[] += h[]*a*d;
	    }
	  }
    }
    boundary ((scalar *){q});
#elif 0
    assert (nl == 1); // fixme
    scalar h = hl[0];
    foreach()
      foreach_dimension()
        if (h[-1] > dry && h[] > dry && h[1] > dry) {
	  double Dp = eta[1] - eta[], Dm = eta[] - eta[-1];
	  if (Dp*Dm < 0. && ((eta[2] + eta[] - 2.*eta[1])*
			     (eta[-1] + eta[1] - 2.*eta[]) < 0. ||
			     (eta[-1] + eta[1] - 2.*eta[])*
			     (eta[-2] + eta[] - 2.*eta[-1]) < 0.)) {
	    double dp, dm;
	    if (fabs(Dp) > fabs(Dm)) {
	      dp = fabs(Dp);
	      dm = fabs(Dm);
	    }
	    else {
	      dp = fabs(Dm);
	      dm = fabs(Dp);
	    }
	    double d = min(dm, dp/2.);
	    double a = Dp > 0. ? 1. : -1.;
	    eta[] += a*d;
	    h[] = max(eta[] - zb[], 0.);
	  }
	}
    boundary ({eta, h});
#endif
  }
  
  /**
  We also use a simple implicit scheme to implement quadratic bottom
  friction i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=10^{-4}$. */
  
  foreach() {
#if ML
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/sq(h[]);
#else
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];    
#endif
    foreach_dimension()
      u.x[] /= a;

    /**
    That is also where we update *hmax*. */

    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
  }
  boundary ({hmax, u});
}

/**
### Snapshots

Every 60 minutes, the $h$, $z_b$ and *hmax* fields are interpolated
bilinearly onto a *n x n* regular grid and written on standard
output. */

event snapshots (t += 15; t <= 400) {
  char name[80];
  sprintf (name, "dump-%g", t);
  dump (name);
}

#if 1
event debugo (t += 5)
  dump();
#endif

event figures (t = 60; t <= 180; t += 60)
{
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  char name[80];
  sprintf (name, "eta-%g.png", t);
  output_ppm (etam, mask = m, min = -1, max = 2, file = name, n = 1024,
	      linear = true, box = {{123,14},{177,55}},
	      opt = "-fill white -opaque black");

  sprintf (name, "level-%g.png", t);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = 5, max = 13, file = name, n = 1024,
	      linear = false, box = {{123,14},{177,55}});

  if (t == 60)
    output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-60.png",
		n = 1024,
		linear = true, box = {{140.6,30.8},{154.6,42}},
		opt = "-fill white -opaque black");
  else if (t == 120)
    output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-120.png",
		n = 1024,
		linear = true, box = {{151.8,27.56},{165.8,38.68}},
		opt = "-fill white -opaque black");
  else if (t == 180)
    output_ppm (etam, mask = m, min = -1, max = 2, file = "zoom-180.png",
		n = 1024,
		linear = true, box = {{158.8,24.25},{172.7,35.3}},
		opt = "-fill white -opaque black");
}

/**
After completion of the simulation, doing

~~~bash
make tsunami/plot.png
~~~

will run gnuplot on these files (using the commands in
[tsunami.plot]()) to produce images such as this one:

![Maximum wave elevation (metres) reached over 10 hours.](tsunami/plot.png)

### Movies

This is done every minute (*t++*). The static variable *fp* is *NULL*
when the simulation starts and is kept between calls (that is what
*static* means). The first time the event is called we set *fp* to a
*ppm2mpeg* pipe. This will convert the stream of PPM images into an
mpeg video using ffmpeg externally. 

We use the *mask* option of *output_ppm()* to mask out the dry
topography. Any part of the image for which *m[]* is negative
(i.e. for which *etam[] < zb[]*) will be masked out. */

event movies (t += 0.5) {
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, mask = m, min = -1, max = 2, 
	      n = 1024, linear = true, file = "eta.mp4");

  /**
  After completion this will give the following animation
  
  ![[Animation](tsunami/eta.mpg) of the wave elevation. Dark blue is -2 metres
  and less. Dark red is +2 metres and more.](tsunami/eta.png)
  
  We also use the *box* option to only output a subset of the domain
  (defined by the lower-left, upper-right coordinates). */
  
  output_ppm (etam, mask = m, min = -5, max = 5, 
              n = 1024, linear = true,
	      box = {{140.4,37.51},{142.44,38.61}},
	      file = "eta-sendai.mp4");

  output_ppm (etam, mask = m, min = -5, max = 5, 
              n = 1024, linear = true,
	      box = {{140.5,36.53},{142.52,37.62}},
	      file = "eta-fukushima.mp4");

  output_ppm (etam, mask = m, min = -5, max = 5, 
              n = 1024, linear = true,
	      box = {{141.32,39.42},{143.44,40.50}},
	      file = "eta-miyako.mp4");

  output_ppm (etam, mask = m, min = -5, max = 5, 
              n = 1024, linear = true,
	      box = {{141.11,38.51},{143.19,39.60}},
	      file = "eta-ofunato.mp4");

  /**
  ![[Animation](tsunami/eta-zoom.mpg) of the wave elevation. Dark blue is 
  -2 metres and less. Dark red is +2 metres and more.](tsunami/eta-zoom.png)
  
  And repeat the operation for the level of refinement...*/

  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = 1024,
	      file = "level.mp4");

  /**
  ![[Animation](tsunami/level.mpg) of the level of refinement. Dark blue is 5
  and dark red is 10.](tsunami/level.png)
  
  ...and for the process id for parallel runs. */
  
#if _OPENMP
  foreach()
    etam[] = pid();
  double tmax = omp_get_max_threads() - 1;
  output_ppm (etam, max = tmax, n = 512, file = "pid.mp4");
#endif // _OPENMP
}

/**
![[Animation](tsunami/pid.mpg) of the OpenMP process id.](tsunami/pid.png)

### Tide gauges

We define a list of file names, locations and descriptions and use the
*output_gauges()* function to output timeseries (for each timestep) of
$\eta$ for each location. */

#include "gauges.h"

/**
As before gnuplot (via [tsunami.plot]()) processes these files to
produce this image:

![Comparison between observed and simulated timeseries (hours) of wave
 elevations (metres) for a selection of tide gauges.](tsunami/gauges.png)

## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event adapt (i++) my_adapt();
