#include "grid/multigrid.h"
#include "poisson.h"
#include "diffusion.h"
#include "run.h"
#include "pnp.h"

#define Volt 1.0
#define DT 0.01

scalar phi[];
scalar Cp[], Cm[];
double dt;
int z[2] = {1,-1};
scalar * sp = {Cp, Cm};

#if 1
const face vector kp[] = {1., 1.};
const face vector km[] = {-1., -1.};
vector * K = {kp, km};
#endif

phi[left]  = dirichlet(Volt);
Cp[left]   = dirichlet (exp(-Volt));
Cm[left]   = dirichlet (exp(Volt));
phi[right] = dirichlet (0.);
Cp[right]  = dirichlet (1.);
Cm[right]  = dirichlet (1.);
 
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

#if 1
  ohmic_flux (sp, z, dt, K);
#else
  ohmic_flux (sp, z, dt); // fixme: this does not work yet
#endif

  for (scalar s in sp)
    diffusion (s, dt);

  /**
  Update the electric potential. */

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

int main() {
  N = 32;
  L0 = 5;
  TOLERANCE = 1e-4;
  run();
}
