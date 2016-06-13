/**
# The 2004 Indian Ocean tsunami

The 2004 Indian Ocean tsunami was caused by a large-scale fault
rupture (> 1000 km) at the Indian–Australian and Eurasian–Andaman
plate boundaries. This example uses the fault model of [Grilli et
al, 2007](/src/references.bib#grilli2007) as initial conditions for a 
Saint-Venant solution of the subsequent tsunami. A similar setup is
discussed in [Popinet, 2011](/src/references.bib#popinet2011).

## Solver setup

The following headers specify that we use spherical coordinates and
the [Saint-Venant solver](/src/saint-venant.h) together with
[(dynamic) terrain reconstruction](/src/terrain.h) and the [Okada
fault model](/src/okada.h). We will use [inputs](/src/input.h) only
when restarting a simulation from a snapshot. */

#include "spherical.h"
#include "saint-venant.h"
#include "terrain.h"
#include "okada.h"
#include "input.h"

/**
We then define a few useful macros and constants. */

int maxlevel = 10;
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HMAXE    5e-2 // error on maximum free surface elevation (5 cm)

/**
The maximum number of levels to use can be set as an argument to the
program. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi(argv[1]);

  /**
  Here we setup the domain geometry. We choose to use metre as length
  unit, so we set the radius of the Earth (required for the [spherical
  coordinates](/src/spherical.h)) in metres. The *x* and *y*
  coordinates are longitude and latitude in degrees, so we set the
  size of the box *L0* and the coordinates of the lower-left corner
  *(X0,Y0)* in degrees. */

  Radius = 6371220.;
  // the domain is 54 degrees squared
  size (54.);
  // centered on 94,8 longitude,latitude
  origin (94. - L0/2., 8. - L0/2.);

  /**
  *G* is the acceleration of gravity required by the Saint-Venant
  solver. This is the only dimensional parameter. We rescale it so that
  time is in minutes. */

  // acceleration of gravity in m/min^2
  G = 9.81*sq(60.);

  /**
  When using a tree (i.e. adaptive) discretisation, we want to start
  with the coarsest grid, otherwise we directly refine to the maximum
  level. Note that *1 << n* is C for $2^n$. */

#if TREE
  // 32^2 grid points to start with
  init_grid (1 << MINLEVEL);
#else // Cartesian
  // 1024^2 grid points
  init_grid (1 << maxlevel);
#endif

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

  run();
}

/**
We declare and allocate another scalar field which will be used to
store the maximum wave elevation reached over time. */

scalar hmax[];

/**
## Boundary conditions

We set the normal velocity component on the left, right and bottom
boundaries to a "radiation condition" with a reference sealevel of
zero. The top boundary is always "dry" in this example so can be left
alone. Note that the sign is important and needs to reflect the
orientation of the boundary. */

u.n[left]   = - radiation(0);
u.n[right]  = + radiation(0);
u.n[bottom] = - radiation(0);

/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two *#if...#else* branches selecting
whether the simulation is being run on an (adaptive) tree or a
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

int adapt() {
#if TREE
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
			    maxlevel, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
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

We then consider two cases, either we restart from an existing
snapshot or we start from scratch. To restart, we could use for example

~~~bash
CFLAGS=-DRESTART make tsunami.tst
~~~

The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */

event init (i = 0)
{
  terrain (zb, "/home/popinet/terrain/etopo2", NULL);
#ifdef RESTART
  input_gfs (file = "snapshot.gfs");
  conserve_elevation();
#else // not RESTART
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
  
  fault (x = 94.57, y = 3.83,
	 depth = 11.4857e3,
	 strike = 323, dip = 12, rake = 90,
	 length = 220e3, width = 130e3,
	 U = 18,
	 iterate = adapt);
#endif // not RESTART
}

/**
The 4 other fault segments are triggered at the appropriate times
(seconds converted to minutes). */

event fault2 (t = 272./60.) {
  fault (x = 93.90, y = 5.22,
	 depth = 11.4857e3,
	 strike = 348, dip = 12, rake = 90,
	 length = 150e3, width = 130e3,
	 U = 23,
	 iterate = adapt);
}

event fault3 (t = 588./60.)
{
  fault (x = 93.21, y = 7.41,
	 depth = 12.525e3,
	 strike = 338, dip = 12, rake = 90,
	 length = 390e3, width = 120e3,
	 U = 12,
	 iterate = adapt);
}

event fault4 (t = 913./60.)
{
  fault (x = 92.60, y = 9.70,
	 depth = 15.12419e3,
	 strike = 356, dip = 12, rake = 90,
	 length = 150e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
}

event fault5 (t = 1273./60.)
{
  fault (x = 92.87, y = 11.70,
	 depth = 15.12419e3,
	 strike = 10, dip = 12, rake = 90,
	 length = 350e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
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
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt speed tn\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g %g %ld\n",
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, perf.speed, grid->tn);

  /**
  We also use a simple implicit scheme to implement quadratic bottom
  friction i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=10^{-4}$. */
  
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
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

event snapshots (t += 60; t <= 600) {

#if !_MPI
  printf ("file: t-%g\n", t);
  output_field ({h, zb, hmax}, stdout, n = 1 << maxlevel, linear = true);
  
  /**
  We also save a snapshot file we can restart from. */

  output_gfs (file = "snapshot.gfs", t = t);
#endif

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
mpeg video. 

We use the *mask* option of *output_ppm()* to mask out the dry
topography. Any part of the image for which *m[]* is negative
(i.e. for which *etam[] < zb[]*) will be masked out. */

event movies (t++) {
  static FILE * fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, mask = m, min = -2, max = 2, n = 512, linear = true);

  /**
  After completion this will give the following animation
  
  ![[Animation](tsunami/eta.mpg) of the wave elevation. Dark blue is -2 metres
  and less. Dark red is +2 metres and more.](tsunami/eta.png)
  
  We also use the *box* option to only output a subset of the domain
  (defined by the lower-left, upper-right coordinates). */
  
  static FILE * fp2 = popen ("ppm2mpeg > eta-zoom.mpg", "w");
  output_ppm (etam, fp2, mask = m, min = -2, max = 2, n = 512, linear = true,
	      box = {{91,6},{100,14}});

  /**
  ![[Animation](tsunami/eta-zoom.mpg) of the wave elevation. Dark blue is 
  -2 metres and less. Dark red is +2 metres and more.](tsunami/eta-zoom.png)
  
  And repeat the operation for the level of refinement...*/

  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = maxlevel, n = 512);

  /**
  ![[Animation](tsunami/level.mpg) of the level of refinement. Dark blue is 5
  and dark red is 10.](tsunami/level.png)
  
  ...and for the process id for parallel runs. */
  
#if _OPENMP
  static FILE * fp3 = popen ("ppm2mpeg > pid.mpg", "w");
  foreach()
    etam[] = tid();
  double tmax = omp_get_max_threads() - 1;
  output_ppm (etam, fp3, max = tmax, n = 512);
#endif // _OPENMP
}

/**
![[Animation](tsunami/pid.mpg) of the OpenMP process id.](tsunami/pid.png)

### Tide gauges

We define a list of file names, locations and descriptions and use the
*output_gauges()* function to output timeseries (for each timestep) of
$\eta$ for each location. */

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

/**
As before gnuplot (via [tsunami.plot]()) processes these files to
produce this image:

![Comparison between observed and simulated timeseries (hours) of wave
 elevations (metres) for a selection of tide gauges.](tsunami/gauges.png) 

### Google Earth KML file

We also generate images and a [Keyhole Markup Language]() file which
can be imported into [Google Earth](https://www.google.com/earth/) to
superpose the evolving wave height field on top of Google Earth
data. */

event kml (t += 15)
{
  static FILE * fp = fopen ("eta.kml", "w");
  if (t == 0)
    fprintf (fp,
	     "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	     "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
	     "  <Folder>\n");
  fprintf (fp,
	   "    <GroundOverlay>\n"
	   "      <TimeSpan>\n"
	   "        <begin>2004-12-26T%02d:%02d:00</begin>\n"
	   "        <end>2004-12-26T%02d:%02d:00</end>\n"
	   "      </TimeSpan>\n"
	   "      <Icon>\n"
	   "	    <href>eta-%g.png</href>\n"
	   "      </Icon>\n"
	   "      <LatLonBox>\n"
	   "	    <north>35</north>\n"
	   "	    <south>-19</south>\n"
	   "	    <east>121</east>\n"
	   "	    <west>67</west>\n"
	   "      </LatLonBox>\n"
	   "    </GroundOverlay>\n",
	   8 + ((int)t)/60, ((int)t)%60,
	   8 + ((int)t + 15)/60, ((int)t + 15)%60,
	   t);
  fflush (fp);
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  char name[80];
  sprintf (name, "eta-%g.png", t);
  output_ppm (etam, file = name, mask = m,
	      min = -2, max = 2, n = 1 << maxlevel, linear = true);
  if (t == 600)
    fprintf (fp, "</Folder></kml>\n");
}

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
