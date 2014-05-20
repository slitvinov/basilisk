/**
# Gouy-Chapman Debye layer */

#include "grid/multigrid.h"
#include "poisson.h"
#include "diffusion.h"
#include "run.h"
#include "pnp.h"

#define Volt 1.0
#define DT 0.01

/**
We assume a fully dissolved binary system with ion $Cp$ and counterion $Cm$
of valence $z$, $|z|=1$).*/

scalar phi[];
scalar Cp[], Cm[];
double dt;
int z[2] = {1,-1};
scalar * sp = {Cp, Cm};

/**
Ions are repelled by the electrode due it positive volume conductivity
while counterions are atracted (negative conductivity), */

#if 1
const face vector kp[] = {1., 1.};
const face vector km[] = {-1., -1.};
vector * K = {kp, km};
#endif

/**
On the left it is placed the charged planar electrode set to a constant
potential $\phi =1$. The concentrations of the positive and negative ions
depend exponentially on the voltage electrode, /*

phi[left]  = dirichlet(Volt);
Cp[left]   = dirichlet (exp(-Volt));
Cm[left]   = dirichlet (exp(Volt));

/**
On the right it is located the liquid bulk. Then, the electrical potential
should be zero. Also, the ion concentrations should match the bulk
concentration i.e */

phi[right] = dirichlet (0.);
Cp[right]  = dirichlet (1.);
Cm[right]  = dirichlet (1.);

/**
Initially, we set the ions concentration to their bulk values imposing
a linear decay of electric potential $phi$ */
 
event init (i = 0)
{
  foreach() {
    phi[] = Volt*(1.-x/5.);
    Cp[] = 1.0; 
    Cm[] = 1.0;
  }
  boundary ({phi, Cp, Cm});
}

event integration (i++) {
  dt = dtnext (t, DT);

  /**
  At each instant, the concentration of each species is updated taking into
  account the ohmic transport. */

#if 1
  ohmic_flux (sp, z, dt, K);
#else
  ohmic_flux (sp, z, dt); // fixme: this does not work yet
#endif

  /**
  Then, the concentrations the thermal diffusion is taken into account, */

  for (scalar s in sp)
    diffusion (s, dt);

  /**
  The electric potential $\phi$ has to be re-calculated since the net
  bulk charge has changed,*/

  scalar rhs[];
  foreach() {
    int i = 0;
    rhs[] = 0.;
    for (scalar s in sp)
      rhs[] -= z[i++]*s[];
  }
  poisson (phi, rhs);
}

event result (t = 3.5) {
  foreach()
    fprintf (stderr, "%g %g %g %g \n", x, phi[], Cp[], Cm[]);
}

/**
## Results

After processing by gnuplot (i.e. using `make debye/plot.png` with
[debye.plot]()) we get

![Profiles of electric potential and concentrations](debye/plot.png) */

int main() {
  N = 32;
  L0 = 5;
  TOLERANCE = 1e-4;
  run();
}
