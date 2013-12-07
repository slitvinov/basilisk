#define GRIDNAME "Cartesian 1D"

#define I     (point.i - 1)
#define J     -0.5
#define DELTA (1./point.n)

struct _Point {
  char * data;
  int i, n;
  int level; // only to return level in locate()
};
#define NN point.n // for output_stencil()

@define data(k,l) ((double *)&point.data[(point.i + k)*datasize + (l) - (l)])
@define allocated(...) true

@define POINT_VARIABLES VARIABLES

@def foreach(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 1; _k <= point.n; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach() } OMP_END_PARALLEL()

@def foreach_boundary(d,corners)
  {
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  {
    point.i = d == right ? point.n : 1;
    POINT_VARIABLES
@
@define end_foreach_boundary() }}

@define foreach_boundary_face(d)    foreach_boundary(d,true)
@define end_foreach_boundary_face() end_foreach_boundary()

@def foreach_boundary_ghost(d) { _OMPSTART /* for face reduction */
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  {
    point.i = (d == right ? point.n : 1) + ig;
    POINT_VARIABLES
@
@define end_foreach_boundary_ghost() } _OMPEND }

void init_grid (int n)
{
  init_solver();
  Point * p = malloc(sizeof(Point));
  size_t len = (n + 2)*datasize;
  p->n = n;
  p->data = malloc (len);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->data;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = p;
  trash (all);
}

void free_grid (void)
{
  delete (all);
  Point * p = grid;
  free (p->data);
  free (p);
  free_solver();
}

void realloc_scalar (void)
{
  Point * p = grid;
  size_t len = (p->n + 2)*datasize;
  p->data = realloc (p->data, len);
  int oldatasize = datasize - sizeof(double);
  char * data = p->data + (p->n + 1)*oldatasize;
  for (int i = p->n + 1; i > 0; i--, data -= oldatasize)
    memmove (data + i*sizeof(double), data, oldatasize);
}

#if TRASH
# undef trash
# define trash cartesian1D_trash
#endif

void cartesian1D_trash (void * alist)
{
  scalar * list = alist;
  Point * p = grid;
  char * data = p->data;
  for (int i = 0; i < p->n + 2; i++, data += datasize) {
    double * v = (double *) data;
    for (scalar s in list)
      v[s] = undefined;
  }
}

#include "cartesian-common.h"

Point locate (double xp, double yp)
{
  Point point = *((Point *)grid);
  point.i = (xp - X0)/L0*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n) ? 0 : - 1;
  return point;
}

void cartesian1D_methods()
{
  cartesian_methods();
}
