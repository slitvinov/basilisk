/* parallel boundary conditions for fluxes */

int main (int argc, char ** argv)
{
  int depth = 4;
  init_grid (2);
  mpi_partitioning();

  origin (-0.6, -0.6);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)), NULL);

  output_cells (stdout);
  
  face vector u[];
  foreach_face() {
    u.x[] = 1.;
    fprintf (stderr, "%g %g\n", x, y);
  }
  boundary_flux ({u});
  fflush (stderr);

  foreach()
    assert (u.x[1,0] - u.x[] + u.y[0,1] - u.y[] == 0.);
}
