#include <math.h>

struct _Data {
  double a, b, res, dp;
};

#include GRID // works with all "multigrid" grids i.e. quadtree.c or multigrid.c
#include "utils.c"
#include "wavelet.c"

void symmetry_level (void * grid, var v, int l)
{
  foreach_boundary_level (grid, right, l)  { val(v,+1,0) = val(v,0,0); } 
  end_foreach_boundary_level();
  foreach_boundary_level (grid, left, l)   { val(v,-1,0) = val(v,0,0); } 
  end_foreach_boundary_level();
  foreach_boundary_level (grid, top, l)    { val(v,0,+1) = val(v,0,0); } 
  end_foreach_boundary_level();
  foreach_boundary_level (grid, bottom, l) { val(v,0,-1) = val(v,0,0); } 
  end_foreach_boundary_level();
}

void relax (void * grid, var a, var b, int l)
{
  foreach_level (grid, l) {
    val(a,0,0) = (val(a,1,0) + val(a,-1,0) +
		  val(a,0,1) + val(a,0,-1) 
		  - delta*delta*val(b,0,0))/4.;
  } end_foreach_level();
  symmetry_level (grid, a, l);
}

void residual (void * grid, var a, var b, var res)
{
  foreach (grid) {
    val(res,0,0) = val(b,0,0) + 
      (4.*val(a,0,0) - val(a,1,0) - val(a,-1,0) - val(a,0,1) - val(a,0,-1))/(delta*delta);
  } end_foreach();
}

void cycle (void * grid, int depth, var a, var res, var dp, int nrelax)
{
  restriction (grid, res);
  foreach_level (grid, 0) { val(dp,0,0) = 0.; } end_foreach_level();
  symmetry_level (grid, dp, 0);
  for (int l = 1; l <= depth; l++) {
    foreach_level (grid, l) {
      /* bilinear interpolation from parent */
      val(dp,0,0) = (9.*coarse(dp,0,0) + 
		     3.*(coarse(dp,childx,0) + coarse(dp,0,childy)) + 
		     coarse(dp,childx,childy))/16.;
    } end_foreach_level();
    symmetry_level (grid, dp, l);
    for (int i = 0; i < nrelax; i++)
      relax (grid, dp, res, l);
  }
  foreach(grid) { val(a,0,0) += val(dp,0,0); } end_foreach();
  symmetry (grid, a);
}

int main()
{
  int depth = 6, nrelax = 4;
  void * grid = init_grid(1 << depth);
  var a = var(a), b = var(b), res = var(res), dp = var(dp);

  foreach(grid) {
    val(b,0,0) = -8.*pi*pi*cos(2.*pi*x)*cos(2.*pi*y);
  } end_foreach();
  
  residual (grid, a, b, res);
  for (int i = 0; i < 20; i++) {
    cycle (grid, depth, a, res, dp, nrelax);
    residual (grid, a, b, res);
    double max = 0.;
    foreach(grid) { if (fabs(val(res,0,0)) > max) max = fabs(val(res,0,0)); } end_foreach();
    fprintf (stderr, "%d %g\n", i, max);
  }

  double max = -1e10, min = 1e10, avg = 0.;
  int n = 0;
  foreach(grid) {
    double e = val(a,0,0) - cos(2.*pi*x)*cos(2.*pi*y);
    if (e > max) max = e; else if (e < min) min = e;
    avg += e; n++;
    printf ("%g %g %g %g %g %g\n", x, y, val(a,0,0), val(b,0,0), val(res,0,0), e);
  } end_foreach();
  fprintf (stderr, "# max error %g\n", max(fabs(max-avg/n),fabs(min-avg/n)));

  free_grid(grid);
}
