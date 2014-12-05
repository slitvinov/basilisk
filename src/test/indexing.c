/* parallel Z-indexing */

int main (int argc, char ** argv)
{
  int depth = 4;
  X0 = Y0 = -0.5;
  init_grid (1);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)), NULL);
  
  scalar reference[];
  int i = 0;
  foreach()
    reference[] = i++;

  mpi_partitioning();

  scalar index[];
  z_indexing (index);

  output_cells (stdout);
  foreach() {
    fprintf (stderr, "%g %g %g %g\n", x, y, index[], reference[]);
    assert (index[] == reference[]);
  }

  free_grid ();
}
