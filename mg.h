/* Multigrid solver */

#include "wavelet.h"

void mg_cycle (void * grid,
	       var a, var res, var dp,
	       void (* relax)    (void * grid, var dp, var res, int depth),
	       void (* boundary) (void * grid, var dp, int depth),
	       int nrelax, int minlevel)
{
  restriction (grid, res, res);
  for (int l = minlevel; l <= depth(grid); l++) {
    if (l == minlevel) {
      foreach_level (grid, l)
	dp(0,0) = 0.;
    }
    else 
      /* bilinear interpolation from coarser level */
      foreach_level (grid, l)
	dp(0,0) = 
	      (9.*coarse(dp,0,0) + 
	       3.*(coarse(dp,childx,0) + coarse(dp,0,childy)) + 
	       coarse(dp,childx,childy))/16.;
    (*boundary) (grid, dp, l);
    for (int i = 0; i < nrelax; i++) {
      (*relax) (grid, dp, res, l);
      (*boundary) (grid, dp, l);
    }
  }
  foreach(grid)
    a(0,0) += dp(0,0);
}
