/**
# Runup of a solitary wave on a conical island

In this test case, we compare experimental data with solutions of the
[Saint-Venant](/src/saint-venant.h) and
[Green-Naghdi](/src/green-naghdi.h) equations. For details on the test
case and experimental data many references are available (e.g. [Hou el
al, 2013](/src/references.bib#hou2013), [Nikolos and Delis,
2009](/src/references.bib#nikolos2009), [Lannes and Marche,
2014](/src/references.bib#lannes2014), [Kazolea et al,
2012](/src/references.bib#kazolea2012) etc...). */

#if SAINT_VENANT
# include "saint-venant.h"
# define MGD 0
#else
# include "green-naghdi.h"
# define MGD mgD.i
#endif

/**
The domain is 27.6 metres squared, resolved with a maximum of 256^2^
grid points. */

#define MAXLEVEL 8
#define MINLEVEL 0

int main()
{
  G = 9.81;
  size (27.6);
  init_grid (1 << MAXLEVEL);
  run();
}

/**
We reproduce case B i.e. a relative amplitude of the solitary wave of
0.09. */

double H0 = 0.32, H = 0.32*0.09, T = 2.45;
#define C sqrt(G*(H0 + H))

double sech2 (double x)
{
  double e = exp (-x);
  double s = 2.*e/(1. + e*e);
  return s*s;
}

double eta0 (double t)
{
  return H*sech2 (sqrt(3.*H/(4.*cube(H0)))*C*(t - T));
}

double uleft (double t)
{
  double e = eta0 (t);
  return C*e/(H0 + e);
}

/**
We use the definition of the solitary wave to impose conditions on the
left side (i.e. the "wave paddle"). */

h[left] = H0 + eta0(t);
u.n[left] = uleft (t);
u.t[left] = 0.;

/**
This is the definition of the topography of the conical island. */

double island (double x, double y)
{
  double r0 = 3.6, r1 = 1.1, H0 = 0.625;
  x -= 25.92/2.; y -= L0/2.;
  double r = sqrt(x*x + y*y);
  return r > r0 ? 0. : r < r1 ? H0 : (r - r0)/(r1 - r0)*H0;
}

/**
On the quadtree we need to make sure that refined cells are
initialised with the correct topography. */

#if QUADTREE
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = island (x, y);
}
#endif

event init (i = 0)
{
#if QUADTREE
  /**
  We setup the refinement function and choose to conserve surface
  elevation rather than volume. */

  zb.refine = refine_zb;
  conserve_elevation();
#endif

  /**
  These are the initial conditions i.e. conical island and
  lake-at-rest. */

  foreach() {
    zb[] = island (x, y);
    h[] = max(0., H0 - zb[]);
  }
}

/**
We define the wave gauges corresponding with the experimental
locations. */

Gauge gauges[] = {
  {"WG3",   6.82, 13.05},
  {"WG6",   9.36, 13.80},
  {"WG9",  10.36, 13.80},
  {"WG16", 12.96, 11.22},
  {"WG22", 15.56, 13.80},
  {NULL}
};

/**
We output snapshots of fields at various times. */

event outputfile (t = {9, 12, 13, 14, 20})
{
  static int nf = 0;
  printf ("file: conical-%d\n", nf);
  output_field ({h,zb}, stdout, N, linear = true);

  /**
  Here we output the level of refinement of the grid. */

  printf ("file: level-%d\n", nf);
  scalar l[];
  foreach()
    l[] = level;
  output_field ({l}, stdout, N);
  nf++;
}

/**
Which gives the following wave (left column) and mesh (right column)
evolutions. 

~~~gnuplot Evolution of free surface (left) and mesh (right)
! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 640,1120 font ",8"
set output 'evolution.png'

unset key
unset xtics
unset ytics
  
dry=1e-4
  
set size ratio -1
set pm3d
set pm3d map
# set contour base
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0,	\
0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )
  
set multiplot layout 4,2 scale 1.35,1.35
splot './conical-0' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-0'
splot './conical-1' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-1'
splot './conical-2' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-2'
splot './conical-3' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-3'
unset multiplot
~~~
*/

event logfile (i++) {

  /**
  We output various diagnostics. */

  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i\n");
  fprintf (stderr, "%g %d %g %g %.8f %g %.4g %g %d\n", 
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD);

  /**
  Here we output surface elevation at the various gauges... */

  output_gauges (gauges, {eta});
}

/**
... and compare it with experimental data for the Green-Naghdi
solver (blue) and the Saint-Venant solver (magenta). The results for
the Green-Naghdi solver are very similar to those of [Kazolea et al,
2012](/src/references.bib#kazolea2012), Figure 14 and [Lannes and
Marche, 2014](/src/references.bib#lannes2014), Figure 20 and
significantly better than the results of the Saint-Venant solver
(see also e.g. Figure 18 of [Hou et al,
2013](/src/references.bib#hou2013) and Figure 19 of [Nikolos and
Delis, 2009](/src/references.bib#nikolos2009)) which have too sharp
fronts.

~~~gnuplot Timeseries of surface elevation. Numerical results and experimental data (symbols).
reset
set term @PNG enhanced size 640,800
set output 'gauges.png'
set multiplot layout 5,1 scale 1.,1.
set xrange [3:20]
t0=22.
h0=0.32
plot '../ts2b.txt' u ($1-t0):($4-0.001) pt 7 ps 0.25 t '',	\
'../conicalsv/WG3' u 1:($2-h0) w l lt 4 t '',			\
'./WG3' u 1:($2-h0) w l lw 2 lt 3 t 'WG3'
plot '../ts2b.txt' u ($1-t0):($6-0.0004) pt 7 ps 0.25 t '',	\
'../conicalsv/WG6' u 1:($2-h0) w l lt 4 t '',			\
'./WG6' u 1:($2-h0) w l lw 2 lt 3 t 'WG6'
plot '../ts2b.txt' u ($1-t0):($7) pt 7 ps 0.25 t '',	\
'../conicalsv/WG9' u 1:($2-h0) w l lt 4 t '',		\
'./WG9' u 1:($2-h0) w l lw 2 lt 3 t 'WG9'
plot '../ts2b.txt' u ($1-t0):($8-0.0017) pt 7 ps 0.25 t '',	\
'../conicalsv/WG16' u 1:($2-h0) w l lt 4 t '',		\
'./WG16' u 1:($2-h0) w l lw 2 lt 3 t 'WG16'
plot '../ts2b.txt' u ($1-t0):($9+0.0015) pt 7 ps 0.25 t '',	\
'../conicalsv/WG22' u 1:($2-h0) w l lt 4 t '',		\
'./WG22' u 1:($2-h0) w l lw 2 lt 3 t 'WG22'
unset multiplot
  
! rm -f conical-?
~~~

We adapt the mesh based on the wavelet-estimated error on the free
surface elevation. */

#if QUADTREE
event adapt (i++) {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  astats s = adapt_wavelet ({eta}, (double[]){3e-4}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif
