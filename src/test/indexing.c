/* parallel Z-indexing */

int main (int argc, char ** argv)
{
  int depth = 4;
  origin (-0.5, -0.5, -0.5);
  init_grid (1);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)),
  	  NULL);
  
  scalar reference[];
  int i = 0;
  foreach()
    reference[] = i++;

  mpi_partitioning();

  scalar index[];
  z_indexing (index);

  output_cells (stdout);
  foreach() {
    fprintf (stderr, "%g %g %g %g %g\n", x, y, z, index[], reference[]);
    assert (index[] == reference[]);
  }
}
