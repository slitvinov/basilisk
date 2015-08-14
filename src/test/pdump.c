// generates the dump file read by restore.c

#include "utils.h"

int main()
{
  int depth = 6;
  origin (-0.5, -0.5, -0.5);
  init_grid (1);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)),
  	  NULL);

#if _MPI  
  mpi_partitioning();
#endif
  
  scalar s[];
  foreach()
    s[] = sin(x)*cos(y);
  boundary ({s});
  restriction ({s});
  
  dump (file = "dump", list = {s});
}
