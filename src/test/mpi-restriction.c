#define BGHOSTS 2

int main (int argc, char * argv[])
{
  origin (-0.6, -0.6, -0.6);
  init_grid (1);
  int depth = argc > 1 ? atoi(argv[1]) : 4;
  refine (level < depth - 2 ||
	  level <= depth*(1. - sqrt(sq(x) + sq(y) + sq(z))),
	  NULL);
  mpi_partitioning();
  
  scalar s[];
  foreach()
    s[] = 1;
  boundary ({s});

  restriction ({s});

  output_cells (stdout);

  // check boundary conditions on leaves
  foreach()
    for (int i = -2; i <= 2; i++)
      for (int j = -2; j <= 2; j++)
	for (int k = -2; k <= 2; k++)
	  assert (s[i,j,k] == 1.);

  // check boundary conditions on levels
  scalar s1[];
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      s1[] = 2;
    boundary_level ({s1}, l);
    foreach_level_or_leaf (l)
      for (int i = -1; i <= 1; i++)
	for (int j = -1; j <= 1; j++)
	  for (int k = -1; k <= 1; k++)
	    assert (s1[i,j,k] == 2.);
  }
  
  // check restriction 
  foreach_fine_to_coarse() {
    fprintf (stderr, "res %g %g %g %g %d\n", x, y, z, s[], level);
    assert (s[] == 1.);
  }

  // check face traversal
  foreach_face() {
    fprintf (stderr, "face %g %g %g %g\n", x, y, z, s[] - s[-1,0]);
    assert (s[] == s[-1]);
  }
}
