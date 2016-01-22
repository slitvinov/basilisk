#include "grid/balance.h"
#include "utils.h"
#include "output.h"
#include "check_restriction.h"

void image()
{
  // const scalar pid[] = pid(); // fixme: this should work
  scalar pid[];
  foreach()
    pid[] = pid();
  output_ppm (pid, n = 128, min = 0, max = npe() - 1);
}

void partition (const char * prog)
{
  // generates reference solution (in ref)
  init_grid (1);
  refine (level < 4, NULL);

  refine (level <= 9 && sq(x - 0.5) + sq(y - 0.5) < sq(0.05), NULL);
  
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

  scalar s[];
  foreach()
    s[] = 1;
  boundary ({s});
  scalar * list = {s};
  
  int refined;
  do {
    refined = 0;
    quadtree->refined.n = 0;
    foreach_leaf()
      if (level <= 9 && sq(x - 0.5) + sq(y - 0.5) < sq(0.05)) {
	refine_cell (point, list, 0, &quadtree->refined);
	refined++;
      }
    mpi_boundary_refine (list);
    mpi_boundary_update();
    balance(0);
    foreach()
      assert (s[] == 1);
    boundary (list);
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    image();
  } while (refined);

  debug_mpi (stderr);

  foreach()
    assert (s[] == 1);
  
  check_restriction (s);
}
