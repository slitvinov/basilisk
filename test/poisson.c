#include <math.h>

struct _Data {
  double a, b, res, dp;
};

#include "quadtree.c"
#include "utils.c"
#include "wavelet.c"

#define foreach_level(grid,l) foreach_cell(grid) { if (level == l || cell.flags & leaf) {
#define end_foreach_level() continue; } } end_foreach_cell()
#define foreach_boundary_level(grid,dir,l) foreach_boundary(grid,dir) {	\
  if (level == l || cell.flags & leaf) {
#define end_foreach_boundary_level() continue; } } end_foreach_boundary()

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
  void * grid = init_grid(64);
  var a = var(a), b = var(b), res = var(res), dp = var(dp);

  foreach(grid) {
    val(b,0,0) = -8.*pi*pi*cos(2.*pi*x)*cos(2.*pi*y);
  } end_foreach();
  
  residual (grid, a, b, res);
  for (int i = 0; i < 20; i++) {
    cycle (grid, 6, a, res, dp, 4);
    residual (grid, a, b, res);
    double max = 0.;
    foreach(grid) { if (fabs(val(res,0,0)) > max) max = fabs(val(res,0,0)); } end_foreach();
    fprintf (stderr, "%d %g\n", i, max);
  }

  double max = 0.;
  foreach(grid) {
    double e = val(a,0,0) - cos(2.*pi*x)*cos(2.*pi*y);
    if (fabs(e) > max) max = fabs(e);
    printf ("%g %g %g %g %g %g\n", x, y, val(a,0,0), val(b,0,0), val(res,0,0), e);
  } end_foreach();
  fprintf (stderr, "# max error %g\n", max);

  free_grid(grid);
}
