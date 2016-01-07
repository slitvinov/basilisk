#include "grid/balance.h"
#include "utils.h"
#include "output.h"

void image()
{
  scalar pid[];
  foreach()
    pid[] = pid();
  output_ppm (pid, n = 128, min = 0, max = npe() - 1);
}

int main (int argc, char * argv[])
{
#if PARTITION // define this to generate reference solution (balance.ref)
  init_grid (1);
  refine (level < 4, NULL);
#else
  init_grid (16);
#endif

  refine (level < 5 && x < 0.315 && y < 0.438 && x > 0.25 && y > 0.375, NULL);
  unrefine (x < 0.25 && y > 0.5, NULL);

  refine (level < 6 && x < 0.315 && y > 0.8, NULL);

#if PARTITION
  mpi_partitioning();
#else
  for (int i = 0; i < 1; i++) {
    image();
    balance();
  }
#endif
  
  image();

  debug_mpi (stderr);
}
