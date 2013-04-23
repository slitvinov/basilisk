#define GRIDNAME "1D Cartesian grid"

#include <stdlib.h>

#define I     (point.i - 1)
#define J     ((point.n - 1)/2.)
#define DELTA (1./point.n)

typedef struct {
  char * data;
  int i, n;
} Point;

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

#define foreach_boundary_level(d,l) foreach_boundary(d)
#define end_foreach_boundary_level() end_foreach_boundary()

#define depth() 0

void init_grid (int n)
{
  Point * p = malloc(sizeof(Point));
  p->n = n;
  p->data = calloc ((n + 2), datasize);
  grid = p;
  init_boundaries (nvar);
  init_events();
}

void free_grid (void)
{
  Point * p = grid;
  free (p->data);
  free (p);
  free_boundaries ();
}

Point locate (double x, double y)
{
  Point point = *((Point *)grid);
  point.i = (x + 0.5)*point.n + GHOSTS;
  return point;
}

#include "cartesian-common.h"
