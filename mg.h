/* Multigrid solver */

#include "wavelet.h"

void mg_cycle (void * grid,
	       scalar a, scalar res, scalar dp,
	       void (* relax)    (void * grid, scalar dp, scalar res, int depth),
	       void (* boundary) (void * grid, scalar dp, int depth),
	       int nrelax, int minlevel)
{
  restriction (grid, res, res);
  for (int l = minlevel; l <= depth(grid); l++) {
    if (l == minlevel) {
      foreach_level (grid, l)
	dp[] = 0.;
    }
    else 
      /* bilinear interpolation from coarser level */
      foreach_level (grid, l)
	dp[] = (9.*coarse(dp,0,0) + 
		3.*(coarse(dp,childx,0) + coarse(dp,0,childy)) + 
		coarse(dp,childx,childy))/16.;
    (*boundary) (grid, dp, l);
    for (int i = 0; i < nrelax; i++) {
      (*relax) (grid, dp, res, l);
      (*boundary) (grid, dp, l);
    }
  }
  foreach(grid)
    a[] += dp[];
}
