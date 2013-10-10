// Hele-Shaw flow solver

#include "advection.h"
#include "poisson.h"

scalar p[], divu[];
vector kmu[];

void coefficients (void); // user-provided
mgstats mgp;              // statistics of the Poisson solver

event pressure (i++, last)
{
  coefficients();
  // pressure field
  mgp = poisson_variable (p, divu, kmu);
  // velocity from pressure gradient
  trash ({u});
  foreach_face()
    u.x[] = kmu.x[]*(p[] - p[-1,0])/Delta;
  boundary_normal ({u});
  boundary_tangent ({u});
}
