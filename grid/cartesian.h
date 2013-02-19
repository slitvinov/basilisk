#define GRIDNAME "Cartesian grid"

#include <stdlib.h>

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

typedef struct {
  char * d;
  int i, j, n;
} Point;

#define data(k,l) ((double *)&point.d[((point.i + k)*(point.n + 2) + (point.j + l))*datasize])

#define foreach(grid) {						\
  Point point = *((Point *)grid);				\
  for (point.i = 1; point.i <= point.n; point.i++)		\
    for (point.j = 1; point.j <= point.n; point.j++) {		\
      VARIABLES
#define end_foreach() }}

#define foreach_boundary(grid,d) {				\
  Point point = *((Point *)grid);				\
  for (int _k = 1; _k <= point.n; _k++) {			\
    point.i = d > left ? _k : d == right ? point.n : 1;		\
    point.j = d < top  ? _k : d == top   ? point.n : 1;		\
    VARIABLES
#define end_foreach_boundary() }}

void * init_grid (int n)
{
  Point * grid = malloc(sizeof(Point));
  grid->n = n;
  grid->d = calloc ((n + 2)*(n + 2), datasize);
  return grid;
}

void free_grid (Point * grid)
{
  free (grid->d);
  free (grid);
}

Point locate (void * grid, double x, double y)
{
  Point point = *((Point *)grid);
  point.i = (x + 0.5)*point.n + GHOSTS;
  point.j = (y + 0.5)*point.n + GHOSTS;
  return point;
}
