#include "fractions.h"

extern scalar * interfaces;
extern staggered vector u;
extern double dt;

foreach_dimension()
static void sweep_x (scalar c, scalar cc)
{
  vector n[];
  scalar alpha[], flux[];
  
  reconstruction (c, n, alpha);
  foreach_face(x) {
    double un = u.x[]*dt/Delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double cf = (c[i,0] <= 0. || c[i,0] >= 1.) ? c[i,0] :
      rectangle_fraction (c[i,0],
			  - s*n.x[i,0], n.y[i,0], alpha[i,0],
			  -0.5, -0.5, s*un - 0.5, 0.5);
    flux[] = u.x[]*cf;
  }
  //    boundary_normal ({flux});
  foreach()
    c[] += dt*(flux[] - flux[1,0] + cc[]*(u.x[1,0] - u.x[]))/Delta;
  boundary ({c});
}

event vof_advection (i = 1; i++)
{
  for (scalar c in interfaces) {
    scalar cc[];
    foreach()
      cc[] = (c[] > 0.5); // eq. 19 from Weymouth & Yue
    // alternate directions
    void (* sweep[2]) (scalar, scalar) = {sweep_x, sweep_y};
    for (int d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc);
  }
}
