/**
# Macros $f_s$ and $f_r$ for the FENE-P model

See [log-conform.h](log-conform.h). */

double L2 = 1.;

static void fenep (double trA, double * nu, double * eta) {
  *eta = 1;
  *nu = 1./(1. - trA/L2);
} 

#define f_r(trA, nu, eta) fenep (trA, &(nu), &(eta))
#define f_s(trA, nu, eta) fenep (trA, &(nu), &(eta))

#include "log-conform.h"

event init (i = 0) {
#if AXI
  double dim = 3;
#else
  double dim = dimension;
#endif  
  scalar trac = trA;
  foreach()
    trac[] = dim*L2/(dim + L2);
}
