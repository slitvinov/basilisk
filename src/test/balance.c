#include "grid/balance.h"
#include "utils.h"
#include "output.h"
#include "check_restriction.h"

void image()
{
  scalar pid[];
  foreach()
    pid[] = pid();
  output_ppm (pid, n = 128, min = 0, max = npe() - 1);
}

void refinement()
{
  refine (level < 5 && x < 0.315 && y < 0.438 && x > 0.25 && y > 0.375, NULL);
  unrefine (x < 0.25 && y > 0.5, NULL);
  refine (level < 6 && x < 0.315 && y > 0.8, NULL);  
}

void partition (const char * prog)
{
  // generates reference solution
  init_grid (1);
  refine (level < 4, NULL);

  refinement();
  
  mpi_partitioning();

  char name[80];
  sprintf (name, "ref-%d", pid());
  FILE * fp = fopen (name, "w");
  debug_mpi (fp);
  fclose (fp);

  MPI_Barrier (MPI_COMM_WORLD);
  if (pid() == 0) {
    sprintf (name, "cat ref-* > ref && cp -f ref %s.ref", prog);
    system (name);
  }
}

int main (int argc, char * argv[])
{
  partition (argv[0]);
  
  init_grid (16);
  refinement();

  scalar s[];
  foreach()
    s[] = 1;

  do image(); while (balance(0));

  debug_mpi (stderr);

  foreach()
    assert (s[] == 1);

  check_restriction (s);
}
