int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  x -= 0.1; y -= 0.1;
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

int main (int argc, char * argv[])
{
  int depth = argc > 1 ? atoi(argv[1]) : 4;
  X0 = Y0 = -0.5;
  init_grid (1);
  while (refine_function (refine_circle, &depth, NULL));
  mpi_partitioning();

  scalar s[];
  foreach()
    s[] = 1;
  boundary ({s});

  output_cells (stdout);

  foreach()
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++) {
	fprintf (stderr, "%g %g %g\n", x + i*Delta, y + j*Delta, s[i,j]);
	assert (s[i,j] == 1.);
      }

  free_grid();
}
