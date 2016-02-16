/**
# Charged column immersed in a dielectric medium

The charges, initially uniformly distributed, tend to accumulate at
the interface due to electrostatic repulsion (charge relaxation).

Since the total charge at the column remains constant, the electric
potential distribution at the outer medium will also remain constant.

If the macro *NAVIER* is set to 1 a true EHD problem is solved
otherwise the problem is purely electrostatic. The ohmic conduction is
calculated implicitly although the conservative explicit scheme could
also be used. */

#define NAVIER 1

#include "ehd_implicit.h"
#if NAVIER
# include "navier-stokes/centered.h"
# include "ehd_stress.h"
#endif
#include "fractions.h"

/**
Far away from the conducting column the electric potential is zero. */

phi[bottom] = dirichlet(0);
phi[top]    = dirichlet(0);
phi[right]  = dirichlet(0);
phi[left]   = dirichlet(0);
#if NAVIER
p[top]      = dirichlet(0);
p[bottom]   = dirichlet(0);
p[left]     = dirichlet(0);
p[right]    = dirichlet(0);
#endif

FILE * fp = NULL;

#define beta 3.
#define cond 3.
#define rhoini 0.5
#define R 0.1
#define circle(x,y) (sq(R) - sq(x) - sq(y))

event init (t = 0) {

  /**
  The conducting media will be defined through the scalar *f*. This will
  be useful also to define the electrical properties of the media and
  the initial distribution of charge density. */

  scalar f[];
  face vector s[];
  vertex scalar psi[];
  foreach_vertex ()
    psi[] = circle(x,y);
  fractions (psi, f, s);
  foreach()
    rhoe[] = rhoini*f[];
  boundary ({rhoe});

  /**
  The electrical conductivity and permittivity are defined on faces.
  The face permittivity is constructed from the center values of the VOF
  scalar *f* by averaging. On the other hand, face conductivity is set
  using the reconstructed surface fractions, */

#if 1
  foreach()
    epsilonc[] = f[]*beta + (1. - f[]);
  boundary ({epsilonc});
#endif
  foreach_face() {
    double T = (f[] + f[-1])/2.;
    epsilon.x[] = T*beta + (1. - T);
    K.x[] = cond*s.x[];
  }
  boundary ((scalar *){epsilon,K});
}

/**
This event checks the conservation of the total charge. */

event chargesum (i += 2) {
  double Q = statsf(rhoe).sum;
  static double Q0;
  fprintf (stderr, "%g %g %g\n", t, Q, i == 0 ? 0 : fabs(Q - Q0));
  if (i == 0)
    Q0 = Q;
  else
    assert (fabs(Q - Q0) < 1e-7);
}

/**
At the final instant, when the charge is fully relaxed on the
interface, the radial electric field and the pressure distribution are
output in order to compare with analytical results. */

event epfield (t = 20) {
  foreach() {
    double Ex = (phi[-1,0] - phi[1,0])/(2*Delta);
    double Ey = (phi[0,-1] - phi[0,1])/(2*Delta);
#if NAVIER
    fprintf (fp, "%g %g %g \n", sqrt(x*x + y*y), sqrt(Ex*Ex + Ey*Ey), p[]);
#else
    fprintf (fp, "%g %g \n", sqrt(x*x + y*y), sqrt(Ex*Ex + Ey*Ey));
#endif
  }
  fflush (fp);
}

int main() {

  /**
  The computational domain spans [-1:1][-1:1]. */

  X0 = Y0 = -1.;
  L0 = 2.;
  DT = 1;
  TOLERANCE = 1e-7;

  /**
  We compute the solution for different levels of refinement. */
  
  for (int LEVEL = 6; LEVEL <= 8; LEVEL++) {
    N = 1 << LEVEL;
    char name[80];
    
    /**
    The name of the file with the results is formed with the
    concatenation of *Er* with the corresponding level. */
    
    sprintf (name, "Er-%d", LEVEL);
    fp = fopen (name, "w");
    run();
    fclose (fp);
  }
}

/**
## Results

~~~gnuplot Outer radial electric field distribution as a function of grid refinement.
set output 'Er.png'
set xlabel 'r'
set ylabel 'Er'
set xrange[0:1]
set yrange[0.:0.03]
set sample 1000
R = 0.1
rhoini = 0.5
E(x) = x < R ? 0 : (R*R*rhoini/2/x)
plot 'Er-6' u 1:2 t "Level = 6",  \
     'Er-7' u 1:2 t "Level = 7",  \
     'Er-8' u 1:2 t "Level = 8",  \
     E(x) w l t "Analytical"
~~~

~~~gnuplot Pressure distribution at instant t=20 for different grid refinement.
set output 'P.png'
set xlabel 'r'
set ylabel 'p'
set xrange[0:0.4]
set yrange[-0.0005:0.00005]
set sample 1000
set key right bottom
R = 0.1
rhoini = 0.5
p(x) = x > R ? 0 : -(R*R*rhoini*rhoini/8)
plot 'Er-6' u 1:3 t "Level = 6",  \
     'Er-7' u 1:3 t "Level = 7",  \
     'Er-8' u 1:3 t "Level = 8",  \
      p(x) w l t "Analytical"
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/cylinder.html)
*/
