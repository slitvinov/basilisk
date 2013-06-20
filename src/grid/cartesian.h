#define GRIDNAME "Cartesian"

#include <stdlib.h>

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

struct _Point {
  char * data;
  int i, j, n;
  int level; // only to return level in locate()
};
#define _n point.n // for output_stencil()

@define data(k,l) ((double *)&point.data[((point.i + k)*(point.n + 2) + \
					  (point.j + l))*datasize])

@define POINT_VARIABLES VARIABLES

@def foreach(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  OMP(omp for schedule(static) clause)
  for (int _k = 1; _k <= point.n; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n; point.j++) {
      POINT_VARIABLES
@
@define end_foreach() }} OMP_END_PARALLEL()

@def foreach_boundary(d,corners)
  OMP_PARALLEL()
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  int _start = 1, _end = point.n;
  /* traverse corners only for top and bottom */
  if (corners && d > left) { _start--; _end++; }
  OMP(omp for schedule(static))
  for (int _k = _start; _k <= _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n : 1;
    point.j = d < top  ? _k : d == top   ? point.n : 1;
    POINT_VARIABLES
@
@define end_foreach_boundary() } OMP_END_PARALLEL()

void init_grid (int n)
{
  init_solver();
  Point * p = malloc(sizeof(Point));
  size_t len = (n + 2)*(n + 2)*datasize;
  p->n = n;
  p->data = malloc (len);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->data;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = p;
}

void realloc_scalar (void)
{
  Point * p = grid;
  size_t oldatasize = datasize - sizeof(double);
  size_t len = (p->n + 2)*(p->n + 2);
  p->data = realloc (p->data, len*datasize);
  char * data = p->data + (len - 1)*oldatasize;
  for (int i = len - 1; i > 0; i--, data -= oldatasize)
    memmove (data + i*sizeof(double), data, oldatasize);  
}

void free_grid (void)
{
  Point * p = grid;
  free (p->data);
  free (p);
  free_solver();
}

Point locate (double xp, double yp)
{
  Point point = *((Point *)grid);
  point.i = (xp - X0)/L0*point.n + 1;
  point.j = (yp - Y0)/L0*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n &&
		 point.j >= 1 && point.j <= point.n) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"
