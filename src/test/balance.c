#include "utils.h"
#include "output.h"
#include "check_restriction.h"
#include "refine_unbalanced.h"

#define BGHOSTS 2

void refinement()
{
  refine_unbalanced (level < 5 && y < 0.315 && (1. - x) < 0.438
  		     && y > 0.25 && (1. - x) > 0.375, NULL);
  unrefine (y < 0.25 && (1. - x) > 0.5, NULL);
  refine_unbalanced (level < 6 && y < 0.315 && (1. - x) > 0.8, NULL);  
}

void partition (const char * prog)
{
  // generates reference solution
  init_grid (1);

  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  quadtree->dirty = true;

  refine_unbalanced (level < 4, NULL);

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

#if 1
  init_grid (16);
#else
  init_grid (1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  quadtree->dirty = true;

  refine_unbalanced (level < 4, NULL);

  mpi_partitioning();
#endif
  
  refinement();
  
  scalar s[];
  face vector u[];
  quadtree_trash ({s,u});
  foreach()
    s[] = 1;
  foreach_face()
    u.x[] = 1;
  boundary ({s, u});
  
  balance(0);
  
  debug_mpi (stderr);

  foreach()
    foreach_neighbor()
      assert (s[] == 1);
  check_restriction (s);

  foreach_face()
    for (int i = -2; i <= 2; i++)
      assert (u.x[0,i] == 1);
}
