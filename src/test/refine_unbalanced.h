#define refine_unbalanced(cond, list) do {				\
  quadtree->refined.n = 0;						\
  int refined;								\
  do {									\
    refined = 0;							\
    foreach_leaf()							\
      if (cond)	{							\
        refine_cell (point, list, 0, &quadtree->refined);		\
	refined++;							\
      }									\
  } while (refined);							\
  mpi_boundary_refine (list);						\
  mpi_boundary_update();						\
} while (0)
