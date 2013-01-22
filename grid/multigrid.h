#define GRIDNAME "Multigrid"

#include <stdio.h>
#include "common.h"

#define GHOSTS 1        // number of ghost layers
#define I (point.i - GHOSTS)
#define J (point.j - GHOSTS)

typedef struct {
  Data ** d;
  int level, depth;
  int i, j, n;
} Point;

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

/***** Multigrid variables and macros *****/
#define aparent(k,l) point.d[point.level-1][((point.i+GHOSTS)/2+k)*(point.n/2+2*GHOSTS) + \
					    (point.j+GHOSTS)/2+l]
#define child(k,l)   point.d[point.level+1][(2*point.i-GHOSTS+k)*2*(point.n + GHOSTS) + \
					    (2*point.j-GHOSTS+l)]
#define MULTIGRID_VARIABLES						\
  int    level = point.level;                                   NOT_UNUSED(level);   \
  int    childx = 2*((point.i+GHOSTS)%2)-1;                     NOT_UNUSED(childx);  \
  int    childy = 2*((point.j+GHOSTS)%2)-1;                     NOT_UNUSED(childy);

/***** Data macros *****/
#define data(k,l)  point.d[point.level][(point.i + k)*(point.n + 2*GHOSTS) + point.j + l]
#define fine(a,k,l) field(child(k,l), a, double)
#define coarse(a,k,l) field(aparent(k,l), a, double)

#define foreach_level(grid,l) {						\
  Point point = *((Point *)grid);					\
  point.level = l; point.n = 1 << point.level;				\
  double delta = 1./point.n; NOT_UNUSED(delta);				\
  for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)		\
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {	\
      MULTIGRID_VARIABLES						\
      VARIABLES								\

#define end_foreach_level() }}

#define foreach(grid) foreach_level(grid,point.depth)
#define end_foreach() end_foreach_level()

#define foreach_boundary_level(grid,d,l) {				\
  Point point = *((Point *)grid);					\
  point.level = l; point.n = 1 << point.level;				\
  double delta = 1./point.n;		NOT_UNUSED(delta);		\
  for (int _k = GHOSTS; _k < point.n + GHOSTS; _k++) {			\
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS; \
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS; \
    MULTIGRID_VARIABLES							\
    VARIABLES

#define end_foreach_boundary_level() }}

#define foreach_boundary(grid,d) foreach_boundary_level(grid,d,point.depth)
#define end_foreach_boundary() }}

#define foreach_fine_to_coarse(grid) {					\
  Point point = *((Point *)grid);					\
  point.level = point.depth - 1; point.n = 1 << point.level;		\
  for (; point.level > 0; point.n /= 2, point.level--) {		\
    double delta = 1./point.n;		NOT_UNUSED(delta);		\
    for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)	\
      for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {	\
        MULTIGRID_VARIABLES						\
        VARIABLES

#define end_foreach_fine_to_coarse() } } }

void * init_grid (int n)
{
  int r = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "multigrid: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    r++;
  }
  Point * m = malloc(sizeof(Point));
  m->depth = r;
  m->d = malloc(sizeof(Point *)*(r + 1));
  for (int l = 0; l <= r; l++)
    m->d[l] = calloc (_size(l), sizeof (Data));
  return m;
}

void free_grid (Point * grid)
{
  for (int l = 0; l <= grid->depth; l++)
    free (grid->d[l]);
  free(grid->d);
}
