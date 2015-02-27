static int coarsen_circle (Point point)
{
  return sq(x - 0.1) + sq(y - 0.1) > sq(0.1);
}

scalar s[];

int main (int argc, char * argv[])
{
  X0 = Y0 = -0.5;
  init_grid (argc > 1 ? atoi(argv[1]) : 32);

  coarsen_function (coarsen_circle, NULL);
  
  output_cells (stdout);
  
  foreach()
    s[] = 1.;
  boundary ({s});

  // check boundary conditions on leaves
  // fixme: should be = -GHOSTS; i <= GHOSTS
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
