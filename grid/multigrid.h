#define GRIDNAME "Multigrid"

#include <stdio.h>
#include <assert.h>

#define I      (point.i - GHOSTS)
#define J      (point.j - GHOSTS)
#define DELTA  (1./point.n)

typedef struct {
  char ** d;
  int level, depth;
  int i, j, n;
} Point;

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
#define data(k,l)  \
  ((double *)&point.d[point.level][((point.i + k)*(point.n + 2*GHOSTS) + \
				    point.j + l)*datasize])

/***** Multigrid variables and macros *****/
#define depth()       (((Point *)grid)->depth)
#define _fine(a,k,l)   \
  ((double *)								\
   &point.d[point.level+1][((2*point.i-GHOSTS+k)*2*(point.n + GHOSTS) + \
			    (2*point.j-GHOSTS+l))*datasize])[a]
#define _coarse(a,k,l) \
  ((double *)								\
   &point.d[point.level-1][(((point.i+GHOSTS)/2+k)*(point.n/2+2*GHOSTS) + \
			    (point.j+GHOSTS)/2+l)*datasize])[a]
#define MULTIGRID_VARIABLES					     \
  int    level = point.level;                   NOT_UNUSED(level);   \
  int    childx = 2*((point.i+GHOSTS)%2)-1;     NOT_UNUSED(childx);  \
  int    childy = 2*((point.j+GHOSTS)%2)-1;     NOT_UNUSED(childy);

#define foreach_level(l,...) 						\
  OMP_PARALLEL()							\
  Point point = *((Point *)grid);					\
  point.level = l; point.n = 1 << point.level;				\
  OMP(omp for schedule(static) __VA_ARGS__)				\
  for (int _k = GHOSTS; _k < point.n + GHOSTS; _k++) {			\
    point.i = _k;							\
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {	\
      MULTIGRID_VARIABLES						\
      VARIABLES								\

#define end_foreach_level() }} OMP_END_PARALLEL()

#define foreach(clause) foreach_level(point.depth, clause)
#define end_foreach()   end_foreach_level()

#define foreach_boundary_level(d,l) 					\
  OMP_PARALLEL()							\
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);		\
  Point point = *((Point *)grid);					\
  point.level = l; point.n = 1 << point.level;				\
  OMP(omp for schedule(static))						\
  for (int _k = GHOSTS; _k < point.n + GHOSTS; _k++) {			\
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS; \
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS; \
    MULTIGRID_VARIABLES							\
    VARIABLES

#define end_foreach_boundary_level() } OMP_END_PARALLEL()

#define foreach_boundary(d)    foreach_boundary_level(d,point.depth)
#define end_foreach_boundary() end_foreach_boundary_level()

#define foreach_fine_to_coarse() {					\
  Point _p = *((Point *)grid);						\
  _p.level = _p.depth - 1; _p.n = 1 << _p.level;			\
  for (; _p.level > 0; _p.n /= 2, _p.level--)				\
    OMP_PARALLEL()							\
    Point point = _p;							\
    OMP(omp for schedule(static))					\
    for (int _k = GHOSTS; _k < point.n + GHOSTS; _k++) {		\
      point.i = _k;							\
      for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {	\
        MULTIGRID_VARIABLES						\
        VARIABLES

#define end_foreach_fine_to_coarse() }} OMP_END_PARALLEL() }

void init_grid (int n)
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
    m->d[l] = calloc (_size(l), datasize);
  grid = m;
  init_boundaries (nvar);
}

void free_grid (void)
{
  Point * m = grid;
  for (int l = 0; l <= m->depth; l++)
    free (m->d[l]);
  free(m->d);
  free_boundaries();
}

Point locate (double x, double y)
{
  Point point = *((Point *)grid);
  point.level = point.depth;
  point.n = 1 << point.level;
  point.i = (x + 0.5)*point.n + GHOSTS;
  point.j = (y + 0.5)*point.n + GHOSTS;
  return point;
}

#include "multigrid-common.h"
