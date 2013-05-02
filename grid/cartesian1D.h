#define GRIDNAME "Cartesian 1D"

#include <stdlib.h>

#define I     (point.i - 1)
#define J     ((point.n - 1)/2.)
#define DELTA (1./point.n)

struct _Point {
  char * data;
  int i, n;
  int level; // only to return level in locate()
};
#define _n point.n // for output_stencil()

#define data(k,l) ((double *)&point.data[(point.i + k)*datasize + (l) - (l)])

#define POINT_VARIABLES VARIABLES

#define foreach(clause)							\
  OMP_PARALLEL()							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  Point point = *((Point *)grid);					\
  OMP(omp for schedule(static) clause)					\
  for (int _k = 1; _k <= point.n; _k++) {				\
    point.i = _k;							\
    POINT_VARIABLES
#define end_foreach() } OMP_END_PARALLEL()

#define foreach_boundary(d) {						\
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);		\
  Point point = *((Point *)grid);					\
  {									\
    point.i = d == right ? point.n : 1;					\
    POINT_VARIABLES
#define end_foreach_boundary() }}

#define foreach_boundary_ghost(d) {					\
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);		\
  Point point = *((Point *)grid);					\
  {									\
    point.i = (d == right ? point.n : 1) + ig;				\
    POINT_VARIABLES
#define end_foreach_boundary_ghost() }}

#define foreach_boundary_level(d,l,...) foreach_boundary(d)
#define end_foreach_boundary_level() end_foreach_boundary()

#define depth() 0

void init_grid (int n)
{
  init_solver();
  Point * p = malloc(sizeof(Point));
  size_t len = (n + 2)*datasize;
  p->n = n;
  p->data = malloc (len);
#if TRASH
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->data;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
#endif
  grid = p;
}

void free_grid (void)
{
  Point * p = grid;
  free (p->data);
  free (p);
  free_solver();
}

Point locate (double x, double y)
{
  Point point = *((Point *)grid);
  point.i = (x + 0.5)*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"

void cartesian1D_methods()
{
  cartesian_methods();
}
