#define GRID "Cartesian grid"

#include <stdlib.h>
#include "utils.h"

#define I (point.i - 1)
#define J (point.j - 1)

typedef struct {
  Data * d;
  int i, j, n;
} Point;

#define data(k,l) point.d[(point.i + k)*(point.n + 2) + (point.j + l)]

#define foreach(grid) {						\
  Point point = *((Point *)grid);				\
  double delta = 1./point.n;		NOT_UNUSED(delta);	\
  for (point.i = 1; point.i <= point.n; point.i++)		\
    for (point.j = 1; point.j <= point.n; point.j++) {		\
      VARIABLES
#define end_foreach() }}

#define foreach_boundary(grid,d) {				\
  Point point = *((Point *)grid);				\
  double delta = 1./point.n;		NOT_UNUSED(delta);	\
  for (int _k = 1; _k <= point.n; _k++) {			\
    point.i = d > left ? _k : d == right ? point.n : 1;		\
    point.j = d < top  ? _k : d == top   ? point.n : 1;		\
    VARIABLES
#define end_foreach_boundary() }}

void * init_grid (int n)
{
  Point * grid = malloc(sizeof(Point));
  grid->n = n;
  grid->d = calloc ((n + 2)*(n + 2), sizeof (Data));
  return grid;
}

void free_grid (Point * grid)
{
  free (grid->d);
  free (grid);
}
