#include "heights.h"

attribute {
  double sigma;
  scalar kappa;
};

event init (i = 0) {
  a = new face vector;
}

event stability (i++) {
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face() {
    if (alpha.x[] > amax) amax = alpha.x[];
    if (alpha.x[] < amin) amin = alpha.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;
  for (scalar c in interfaces)
    if (c.sigma) {
      double dtmax = sqrt (rhom*cube(dmin)/(pi*c.sigma));
      if (dtmax < DT)
	DT = dtmax;
    }
}

event acceleration (i++)
{
  scalar * list = NULL;
  for (scalar c in interfaces)
    if (c.sigma) {
      list = list_add (list, c);
      if (!c.kappa)
	c.kappa = new scalar;
      curvature (c, c.kappa);
    }

  foreach_face() {
    a.x[] = 0.;
    for (scalar c in list)
      if (c[] != c[-1,0]) {
	scalar kappa = c.kappa;
	a.x[] = alpha.x[]*c.sigma*(kappa[] + kappa[-1,0])/2.*
	  (c[] - c[-1,0])/Delta;
      }
  }

  free (list);
}
