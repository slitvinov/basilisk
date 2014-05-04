#include "utils.h"

int main ()
{
  init_grid (64);
#if _MPI
  mpi_partitioning();
#endif

  scalar s[];
  foreach()
    s[] = x + y;
  boundary ({s});

  // statsf() uses MPI reduction operations
  stats stat = statsf (s);
  fprintf (stderr, "%g %g %g\n", stat.min, stat.sum, stat.max);
}
