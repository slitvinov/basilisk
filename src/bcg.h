void fluxes_upwind_bcg (const scalar f, const vector u, 
			vector flux,
			double dt)
{
  vector g[];
  gradients ({f}, {g});
  
  trash ({flux});
  foreach_face() {
    double un = dt*u.x[]/delta, s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i,0] + s*min(1., 1. - s*un)*g.x[i,0]*delta/2.;
    double vn = u.y[i,0] + u.y[i,1];
    double fyy = vn < 0. ? f[i,1] - f[i,0] : f[i,0] - f[i,-1];
    f2 -= dt*vn*fyy/(4.*delta);
    flux.x[] = f2*u.x[];
  }
  boundary_normal ({flux});
}
