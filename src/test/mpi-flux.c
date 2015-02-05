/* parallel boundary conditions for fluxes */

int main (int argc, char ** argv)
{
  int minlevel = argc > 1 ? atoi(argv[1]) : 5;

  init_grid (1);
  refine (level <= minlevel*(1. - sqrt(sq((x - 0.5) - 0.1) +
				       sq((y - 0.5)- 0.1))), NULL);

  mpi_partitioning();

  output_cells (stdout);

  scalar s[];
  foreach()
    s[] = noise();
  boundary ({s});
  
  face vector g[];
  foreach_face()
    g.x[] = (s[] - s[-1,0])/Delta;
  boundary_flux ({g});

  foreach_face(x)
    fprintf (stderr, "%g %g %g\n", x, y, g.x[]);
  
  double sum = 0.;
  foreach (reduction(+:sum))
    sum += Delta*(g.x[1,0] - g.x[] + g.y[0,1] - g.y[]);
  fprintf (stderr, "sum: %g\n", sum);
  assert (fabs(sum) < 1e-10);
}
