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
#if _MPI
  mpi_partitioning();
#endif

  scalar s[];
  foreach()
    s[] = 1;
  boundary ({s});

  restriction ({s});

  output_cells (stdout);

  // check boundary conditions on leaves
  foreach()
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	assert (s[i,j] == 1.);

  // check boundary conditions on levels
  scalar s1[];
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      s1[] = 2;
    boundary_level ({s1}, l);
    foreach_level_or_leaf (l)
      for (int i = -1; i <= 1; i++)
	for (int j = -1; j <= 1; j++)
	  assert (s1[i,j] == 2.);
  }
    
  // check restriction 
  foreach_fine_to_coarse() {
    fprintf (stderr, "%g %g %g %d\n", x, y, s[], level);
    assert (s[] == 1.);
  }

  free_grid();
}
