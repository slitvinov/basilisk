int main (int argc, char * argv[])
{
  X0 = Y0 = -0.5;
  init_grid (2);
  mpi_partitioning();

  int depth = argc > 1 ? atoi(argv[1]) : 6;
  refine((level <= depth && x <= -0.25 && y < 0 && y >= -0.25) ||
	 (level <= depth - 1 && y < 0), NULL);
  
  output_cells (stdout);
  
  scalar s[];
  foreach()
    s[] = 1.;
  boundary ({s});

  // check boundary conditions on leaves
  foreach()
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	assert (s[i,j] == 1.);

  // rebalancing
  int nf = 0;
  foreach(reduction(+:nf))
    nf++;
  int npe;
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  nf = max(1, nf/npe);
  scalar index[];
  z_indexing (index);
  foreach()
    fprintf (stderr, "%g %g %d\n", x, y, min(npe - 1, (int)(index[]/nf)));
}