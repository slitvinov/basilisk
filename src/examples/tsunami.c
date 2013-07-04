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

#include "saint-venant.h"
#include "terrain.h"
#include "okada.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 10
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HMAXE    5e-2 // error on maximum free surface elevation (5 cm)

// metres to degrees
double mtd = 360./40075e3;

/**
The following function will be called by the Saint-Venant solver
to setup the initial simulation grid. When using a quadtree
(i.e. adaptive) discretisation, we want to start with the coarsest
grid, otherwise we directly refine to the maximum level. Note that `1
<< n` is C for $2^n$. */

void parameters()
{
#if QUADTREE
  // 32^2 grid points to start with
  N = 1 << MINLEVEL;
#else // Cartesian
  // 1024^2 grid points
  N = 1 << MAXLEVEL;
#endif

/**
Here we setup the domain geometry. For the moment Basilisk only
supports square domains. For this example we have to use degrees as
horizontal units because that is what the topographic database uses
(eventually coordinate mappings will give more flexibility). We set
the size of the box `L0` and the coordinates of the lower-left corner
`(X0,Y0)`. */

  // the domain is 54 degrees squared
  L0 = 54.;
  // centered on 94,8 longitude,latitude
  X0 = 94 - L0/2.;
  Y0 = 8. - L0/2.;

/**
`G` is the acceleration of gravity required by the Saint-Venant
solver. This is the only dimensional parameter. We rescale it so that
time is in minutes, horizontal distances in degrees and vertical
distances in metres. This is a trick to circumvent the current lack of
coordinate mapping. */

  // acceleration of gravity in degrees^2/min^2/m
  G = 9.81*sq(mtd)*sq(60.);
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

u.x[left]   = - radiation(0);
u.x[right]  = + radiation(0);
u.y[bottom] = - radiation(0);

/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two `#if...#else` branches selecting
whether the simulation is being run on an (adaptive) quadtree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height `hmax` -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
`η`) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

/**
We can now use wavelet adaptation on the list of scalars `{η,hmax}`
with thresholds `{ETAE,HMAXE}`. The compiler is not clever enough yet
and needs to be told explicitly that this is a list of `double`s,
hence the `(double[])`
[type casting](http://en.wikipedia.org/wiki/Type_conversion). 

The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({eta, hmax}, (double[]){ETAE,HMAXE},
			    MAXLEVEL, MINLEVEL);
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
[`xyz2kdt` manual](http://gfs.sourceforge.net/wiki/index.php/Xyz2kdt)
for explanations on how to do this.

The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */

void init()
{
  terrain (zb, "/home/popinet/terrain/etopo2", NULL);
  conserve_elevation();

/**
The initial still water surface is at $z=0$ so that the water depth $h$ is... */

  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});

/**
The initial deformation is given by an Okada fault model with the
following parameters. The `iterate = adapt` option will iterate this
initialisation until our `adapt()` function above returns zero
i.e. until the deformations are resolved properly. */

  fault (x = 94.57, y = 3.83,
	 depth = 11.4857e3,
	 strike = 323, dip = 12, rake = 90,
	 length = 220e3, width = 130e3,
	 U = 18,
	 iterate = adapt);
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

We output simple summary statistics for `h` and `u.x` on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);

/**
We also use a simple implicit scheme to implement quadratic bottom
friction i.e.
$$
\frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
$$
with $C_f=10^{-4}$. */

  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/(h[]*mtd);
    foreach_dimension()
      u.x[] /= a;

/**
That is also where we update `hmax`. */

    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
  }
  boundary ({hmax, u});
}

/**
### Snapshots

Every 60 minutes, the $h$, $z_b$ and `hmax` fields are interpolated
bilinearly onto a `n x n` regular grid and written on standard
output. */

event snapshots (t += 60; t <= 600) {
  printf ("file: t-%g\n", t);
  output_field ({h, zb, hmax}, stdout, n = 1 << MAXLEVEL, linear = true);
}

/**
After completion of the simulation, doing

~~~bash
make tsunami.png
~~~

will run gnuplot on these files (using the commands in
[tsunami.plot]()) to produce images such as this one:

![Maximum wave elevation (metres) reached over 10 hours.](tsunami.png)

### Movies

This is done every minute (`t++`). The static variable `fp` is `NULL`
when the simulation starts and is kept between calls (that is what
`static` means). The first time the event is called we set `fp` to a
`ppm2mpeg` pipe. This will convert the stream of PPM images into an
mpeg video using ffmpeg externally. 

We use the `mask` option of `output_ppm()` to mask out the dry
topography. Any part of the image for which `m[]` is negative
(i.e. for which `etam[] < zb[]`) will be masked out. */

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

/**
After completion this will give the following animation

![[Animation](eta.mpg) of the wave elevation. Dark blue is -2 metres
 and less. Dark red is +2 metres and more.](eta.png)

We also use the `box` option to only output a subset of the domain
(defined by the lower-left, upper-right coordinates). */

  static FILE * fp2 = NULL;
  if (!fp2) fp2 = popen ("ppm2mpeg > eta-zoom.mpg", "w");
  output_ppm (etam, fp2, mask = m, min = -2, max = 2, n = 512, linear = true,
	      box = {{89,8},{98,16}});

/**
![[Animation](eta-zoom.mpg) of the wave elevation. Dark blue is -2 metres
 and less. Dark red is +2 metres and more.](eta-zoom.png)

And repeat the operation for the level of refinement...*/

  static FILE * fp1 = NULL;
  if (!fp1) fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = MAXLEVEL, n = 512);

/**
![[Animation](level.mpg) of the level of refinement. Dark blue is 5
 and dark red is 10.](level.png)

...and for the process id for parallel runs. */

#if _OPENMP
  static FILE * fp3 = NULL;
  if (!fp3) fp3 = popen ("ppm2mpeg > pid.mpg", "w");
  foreach()
    etam[] = pid();
  double tmax = omp_get_max_threads() - 1;
  output_ppm (etam, fp3, max = tmax, n = 512);
#endif // _OPENMP
}

/**
![[Animation](pid.mpg) of the OpenMP process id.](pid.png)

### Tide gauges

We define a list of file names, locations and descriptions and use the
`output_gauges()` function to output timeseries (for each timestep) of
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
 elevations (metres) for a selection of tide gauges.](tsunami_gauges.png)

## Adaptivity

We apply our `adapt()` function at every timestep. */

event do_adapt (i++) adapt();

/**
And finally we call the `run()` method of the Saint-Venant solver to
perform the integration. */

int main() { run(); }
