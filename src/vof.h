#include "fractions.h"

extern scalar * interfaces;
extern staggered vector uf;
extern double dt;

event defaults (i = 0)
{
  CFL = 0.5;
#if QUADTREE
  for (scalar c in interfaces) {
    c.prolongation = fraction_prolongation;
    c.refine = fraction_refine;
  }
#endif
}

foreach_dimension()
static void sweep_x (scalar c, scalar cc)
{
  vector n[];
  scalar alpha[], flux[];
  
  reconstruction (c, n, alpha);
  foreach_face(x) {
    double un = uf.x[]*dt/Delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double cf = (c[i,0] <= 0. || c[i,0] >= 1.) ? c[i,0] :
      rectangle_fraction (- s*n.x[i,0], n.y[i,0], alpha[i,0],
			  -0.5, -0.5, s*un - 0.5, 0.5);
    flux[] = uf.x[]*cf;
  }

#if QUADTREE
  // like boundary_normal() but for a single dimension
  foreach_halo_fine_to_coarse() {
    flux[] = (fine(flux,0,0) + fine(flux,0,1))/2.;
    if (is_leaf (neighbor(1,0)))
      flux[1,0] = (fine(flux,2,0) + fine(flux,2,1))/2.;
  }
#endif

  foreach()
    c[] += dt*(flux[] - flux[1,0] + cc[]*(uf.x[1,0] - uf.x[]))/Delta;
  boundary ({c});
}

event vof (i = 1; i++)
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
