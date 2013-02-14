void parameters ()
{
  // number of grid points
  N = 64;
  // viscosity
  NU = 1e-3;
  // end time
  TMAX = 300;
  // maximum timestep
  DT = 1e-1;
  // CFL number
  CFL = 0.8;
}

void initial_conditions (void * grid)
{
  /* default to zero */
}

void boundary_p (void * grid, var p, int l)
{
  foreach_boundary_level (grid, right, l)  p(1,0)  = p(0,0);
  foreach_boundary_level (grid, left, l)   p(-1,0) = p(0,0);
  foreach_boundary_level (grid, top, l)    p(0,1)  = p(0,0);
  foreach_boundary_level (grid, bottom, l) p(0,-1) = p(0,0);
}

void boundary_u (void * grid, var u, var v)
{
  foreach_boundary (grid, right) {
    u(1,0) = 0.;
    v(1,0) = - v(0,0);
  }
  foreach_boundary (grid, left) {
    u(-1,0) = - u(1,0);
    u(0,0) = 0.;
    v(-1,0) = - v(0,0);
  }
  foreach_boundary (grid, top) {
    v(0,1) = 0.;
    u(0,1) = 2. - u(0,0);
  }
  foreach_boundary (grid, bottom) {
    v(0,-1) = - v(0,1);
    v(0,0) = 0.;
    u(0,-1) = - u(0,0);
  }
}

double energy (void * grid)
{
  double se = 0.;
  foreach (grid)
    se += (sq(u(0,0) + u(1,0)) + sq(v(0,0) + v(0,1)))/8.*delta*delta;
  return se*L0*L0;
}

void output_field (void * grid, var f, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  double delta = 1./N;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double x = delta*i - 0.5 + delta/2., y = delta*j - 0.5 + delta/2.;
      fprintf (fp, "%g %g %g\n", x, y, interpolate (grid, f, x, y));
    }
    fputc ('\n', fp);
  }
}

var un = var(un); /* we need another variable */

int events (void * grid, int i, double t, double dt)
{
  if (i % 10 == 0) {
    double du = change (grid, u, un);
    if (i > 0 && du < 1e-4)
      return 1; /* stop */
    fprintf (stderr, "t: %f %.9f %g\n", t, energy (grid), du);
  }
  if (i % 100 == 0)
    output_field (grid, u, stdout);
  return 0; /* continue */
}

void end (void * grid)
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (grid, u, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (grid, v, x, 0));
  fclose (fp);
}
