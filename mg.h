/* Multigrid solver */

void mg_cycle (scalar a, scalar res, scalar dp,
	       void (* relax)    (scalar dp, scalar res, int depth),
	       void (* boundary) (scalar dp, int depth),
	       int nrelax, int minlevel)
{
  restriction (res, res);
  for (int l = minlevel; l <= depth(); l++) {
    if (l == minlevel) {
      foreach_level (l)
	dp[] = 0.;
    }
    else 
      /* bilinear interpolation from coarser level */
      foreach_level (l)
	dp[] = (9.*coarse(dp,0,0) + 
		3.*(coarse(dp,childx,0) + coarse(dp,0,childy)) + 
		coarse(dp,childx,childy))/16.;
    (*boundary) (dp, l);
    for (int i = 0; i < nrelax; i++) {
      (*relax) (dp, res, l);
      (*boundary) (dp, l);
    }
  }
  foreach()
    a[] += dp[];
}
