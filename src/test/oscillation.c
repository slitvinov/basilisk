/**
# Shape oscillation of an inviscid droplet

This test case is discussed in [Popinet, 2009](/src/references.bib#popinet2009).

A two-dimensional elliptical droplet (density ratio 1/1000) is
released in a large domain. Under the effect of surface-tension the
shape of the droplet oscillates around its (circular) equilibrium
shape. The fluids inside and outside the droplet are inviscid so
ideally no damping of the oscillations should occur. As illustrated on
the figures some damping occurs in the simulation due to numerical
dissipation.

This simulation is also a stringent test case of the accuracy of the
surface tension representation as no explicit viscosity can damp
eventual parasitic currents. 

We use either the momentum-conserving or standard Navier--Stokes
solver with VOF interface tracking and surface tension. */

#if MOMENTUM

# include "momentum.h"
# include "tension.h"
# define cf

#else // standard centered Navier--Stokes solver

# include "navier-stokes/centered.h"
# include "vof.h"
# include "tension.h"

/**
The interface is represented by the volume fraction field *f*. */

scalar f[], * interfaces = {f};

/**
The density inside the droplet is one and outside 10^-3^. */

#define rho(f) (clamp(f,0,1)*(1. - 1e-3) + 1e-3)

/**
We have the option of using some "smearing" of the density jump. */

#if 0
#define cf f
#else
scalar cf[];
#endif

/**
The density is variable. We allocate a new field to store its
inverse. */

face vector alphav[];

/**
The density is defined at each timestep via the *properties()* event
declared by the Navier--Stokes solver. */

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *cf* with the
  vertex-average of *f*. */

#ifndef cf
  foreach()
    cf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
  boundary ({cf});
#endif

  /**
  The inverse of the density $\alpha$ is then given by the
  face-averaged value of *cf* and the arithmetic average of density
  defined by *rho()*. */

  foreach_face() {
    double cm = (cf[] + cf[-1])/2.;
    alphav.x[] = 1./rho(cm);
  }
  boundary ((scalar *){alphav});  
}

#endif // standard centered Navier--Stokes solver

/**
The diameter of the droplet is 0.2. */

#define D 0.2

/**
We will vary the level of refinement to study convergence. */

FILE * fp = NULL;
int LEVEL;

int main() {

  /**
  The density is variable. */

#if MOMENTUM
  rho1 = 1, rho2 = 1e-3;
#else
  alpha = alphav;
#endif
  
  /**
  The surface tension is unity. Decreasing the tolerance on the
  Poisson solve improves the results. We cleanup existing files and
  vary the level of refinement. */

  f.sigma = 1.;
  TOLERANCE = 1e-4;
  remove ("error");
  remove ("laplace");
  for (LEVEL = 5; LEVEL <= 8; LEVEL++) {
    N = 1 << LEVEL;
    
    /**
    We open a file indexed by the level to store the time evolution of
    the kinetic energy. */

    char name[80];
    sprintf (name, "k-%d", LEVEL);
    fp = fopen (name, "w");
    run();
    fclose (fp);
  }

  /**
  We use *grep* to filter the lines generated by gnuplot containing
  the results of the fits (see below). */

  system ("grep ^fit out >> log");
}

event init (i = 0) {

  /**
  We initialise the shape of the interface, a slightly elliptic droplet. */

  fraction (f, D/2.*(1. + 0.05*cos(2.*atan2(y,x))) - sqrt(sq(x) + sq(y)));

#ifndef cf
  foreach()
    cf[] = f[];
#endif
}

/**
At each timestep we output the kinetic energy. */

event logfile (i++; t <= 1) {
  double ke = 0.;
  foreach (reduction(+:ke))
#if MOMENTUM
    ke += sq(Delta)*(sq(q.x[]) + sq(q.y[]))/rho[];
#else
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(cf[]);
#endif
  fprintf (fp, "%g %g %d\n", t, ke, mgp.i);
  fflush (fp);
}

/**
At the end of the simulation, we use gnuplot to fit a function of the form
$$
k(t) = ae^{-bt}(1-\cos(ct))
$$
to the kinetic energy. This gives estimates of the oscillation
pulsation *c* and of the damping *b*.

We also compute the relative error on the pulsation, using the
theoretical value $\omega_0$ as reference. */

event fit (t = end) {
  FILE * fp = popen ("gnuplot 2>&1", "w");
  fprintf (fp, 
           "k(t)=a*exp(-b*t)*(1.-cos(c*t))\n"
           "a = 3e-4\n"
           "b = 1.5\n"
	   "\n"
           "D = %g\n"
           "n = 2.\n"
           "sigma = 1.\n"
           "rhol = 1.\n"
           "rhog = 1./1000.\n"
           "r0 = D/2.\n"
           "omega0 = sqrt((n**3-n)*sigma/((rhol+rhog)*r0**3))\n"
	   "\n"
           "c = 2.*omega0\n"
           "fit k(x) 'k-%d' via a,b,c\n"
	   "level = %d\n"
	   "res = D/%g*2.**level\n"
           "print \"fit \", res, a, b, c, D\n"
	   "\n"
	   "set table 'fit-%d'\n"
	   "plot [0:1] 2.*a*exp(-b*x)\n"
	   "unset table\n"
	   "\n"
	   "set print 'error' append\n"
	   "print res, c/2./omega0-1., D\n"
	   "\n"
	   "set print 'laplace' append\n"
	   "empirical_constant = 30.\n"
	   "print res, (1./(b**2.*D**3.))*empirical_constant**2, D\n"
	   "\n",
	   D, LEVEL, LEVEL, L0, LEVEL);
  pclose (fp);
}

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D oscillation.gfv", "w");
  output_gfs (fp);
}
#endif

#if TREE
event adapt (i++) {
#if MOMENTUM
  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary ((scalar *){u});
#endif
  adapt_wavelet ({f,u}, (double[]){5e-3,1e-3,1e-3}, LEVEL);
}
#endif

/**
## Results

~~~gnuplot Evolution of the kinetic energy as a function of time for the spatial resolutions (number of grid points per diameter) indicated in the legend. The black lines are fitted decreasing exponential functions.
set xlabel 'Time'
set ylabel 'Kinetic energy'
set logscale y
plot [0:1][8e-5:]'k-8' t "51.2" w l, 'k-7' t "25.6" w l,               \
  'k-6' t "12.8" w l, 'k-5' t "6.4" w l,			       \
  'fit-8' t "" w l lt 7, 'fit-7' t "" w l lt 7, 'fit-6' t "" w l lt 7, \
  'fit-5' t "" w l lt 7
~~~

~~~gnuplot Relative error in the oscillation frequency as a function of resolution.
set xlabel 'Diameter (grid points)'
set ylabel 'Frequency error (%)'
set logscale x 2
unset grid
set xzeroaxis
set key spacing 1.5 top right
ftitle(a,b,c) = sprintf("%.0f/x^{%4.2f} (%s)", exp(a), -b, c)
f(x)=a+b*x
fit f(x) 'error' u (log($1)):(log(abs($2)*100.)) via a,b
f1(x)=a1+b1*x
fit f1(x) '../oscillation-momentum/error' u (log($1)):(log(abs($2)*100.)) via a1,b1
plot 'error' u ($1):(abs($2)*100.) t "" w p pt 5 ps 2, \
     '../oscillation-momentum/error' u ($1):(abs($2)*100.) t "" w p pt 7 ps 2, \
      exp(f(log(x))) t ftitle(a,b,"standard"), \
      exp(f1(log(x))) t ftitle(a1,b1,"momentum")
~~~

The amount of numerical damping can be estimated by computing an
equivalent viscosity. With viscosity, kinetic energy is expected to
decrease as:
$$\exp(-C\nu/D^2t)$$
where $C$ is a constant, $\nu$ the viscosity and $D$ the droplet
diameter. Using curve fitting the damping coefficient $b=C\nu/D^2$
can be estimated (black curves on Figure \ref{kinetic}). An
equivalent Laplace number can then be computed as:
$$La=\frac{\sigma D}{\rho\nu^2}=\frac{\sigma C^2}{\rho b^2 D^3}$$
The equivalent Laplace number depends on spatial resolution as
illustrated below.

~~~gnuplot Equivalent Laplace number estimated from the numerical damping of kinetic energy.
set xlabel 'Diameter (grid points)'
set ylabel 'Equivalent Laplace number'
set grid
set key bottom right
plot 'laplace' t "standard" w p pt 5 ps 2, \
     '../oscillation-momentum/laplace' t "momentum" w p pt 7 ps 2
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/oscillation.html)
*/
