/* parallel Z-indexing */

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  x -= 0.1; y -= 0.1;
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

int main (int argc, char ** argv)
{
  int depth = 4;
  X0 = Y0 = -0.5;
  init_grid (1);
  while (refine_function (refine_circle, &depth, NULL));
  
  scalar reference[];
  int i = 0;
  foreach()
    reference[] = i++;

#if _MPI
  mpi_partitioning();
#endif

  scalar index[];
  z_indexing (index);

  output_cells (stdout);
  foreach() {
    fprintf (stderr, "%g %g %g %g\n", x, y, index[], reference[]);
    assert (index[] == reference[]);
  }

  free_grid ();
}
