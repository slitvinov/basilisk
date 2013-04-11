#define GRIDNAME "Cartesian grid"

#include <stdlib.h>

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

typedef struct {
  char * data;
  int i, j, n;
} Point;

#define data(k,l) ((double *)&point.data[((point.i + k)*(point.n + 2) + \
					  (point.j + l))*datasize])

#define foreach(clause)							\
  OMP_PARALLEL()							\
  Point point = *((Point *)grid);					\
  OMP(omp for schedule(static) clause)					\
  for (int _k = 1; _k <= point.n; _k++) {				\
    point.i = _k;							\
    for (point.j = 1; point.j <= point.n; point.j++) {			\
      VARIABLES
#define end_foreach() }} OMP_END_PARALLEL()

#define foreach_boundary(d) 						\
  OMP_PARALLEL()							\
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);		\
  Point point = *((Point *)grid);					\
  OMP(omp for schedule(static))						\
  for (int _k = 1; _k <= point.n; _k++) {				\
    point.i = d > left ? _k : d == right ? point.n : 1;			\
    point.j = d < top  ? _k : d == top   ? point.n : 1;			\
    VARIABLES
#define end_foreach_boundary() } OMP_END_PARALLEL()

#define foreach_boundary_level(d,l) foreach_boundary(d)
#define end_foreach_boundary_level() end_foreach_boundary()

#define depth() 0

void init_grid (int n)
{
  Point * p = malloc(sizeof(Point));
  p->n = n;
  p->data = calloc ((n + 2)*(n + 2), datasize);
  grid = p;
  init_boundaries (nvar);
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
  point.j = (y + 0.5)*point.n + GHOSTS;
  return point;
}
