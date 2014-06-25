#define GRIDNAME "Multigrid 1D"

#define I      (point.i - GHOSTS)
#define J      -0.5
#define DELTA  (1./point.n)

struct _Point {
  char ** d;
  int i, j, n;
  int level, depth;
};
#define NN point.n // for output_stencil()

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
@def data(k,l)
((double *)&point.d[point.level][(point.i + k)*datasize + (l) - (l)]) @
@define allocated(...) true

/***** Multigrid variables and macros *****/
@define depth()       (((Point *)grid)->depth)
@def fine(a,k,l)
  ((double *)
   &point.d[point.level+1][(2*point.i-GHOSTS+k)*datasize + (l) - (l)])[a]
@
@def coarse(a,k,l)
  ((double *)
   &point.d[point.level-1][((point.i+GHOSTS)/2+k)*datasize + (l) - (l)])[a]
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 0
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
@
@def foreach_level(l)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = l; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach_level() } OMP_END_PARALLEL()

@def foreach(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = point.depth; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach() } OMP_END_PARALLEL()

@def foreach_boundary_level(d,l,corners) {
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = l; point.n = 1 << point.level;
  {
    point.i = d == right ? point.n : 1;
    POINT_VARIABLES
@
@define end_foreach_boundary_level() }}

@define foreach_boundary_face(d) foreach_boundary_level(d,point.depth,true)
@define end_foreach_boundary_face() end_foreach_boundary_level()

@def foreach_boundary(d,corners) foreach_boundary_level(d,point.depth,corners) @
@def end_foreach_boundary() end_foreach_boundary_level() @

@def foreach_boundary_ghost(d) { _OMPSTART /* for face reduction */
  int ig = _ig[d], jg = _jg[d];	NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = point.depth; point.n = 1 << point.level;
  {
    point.i = (d == right ? point.n : 1) + ig;
    POINT_VARIABLES
@
@define end_foreach_boundary_ghost() } _OMPEND }

@define is_coarse() (point.level < depth())

@def foreach_fine_to_coarse() {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point _p = *((Point *)grid);
  _p.level = _p.depth - 1; _p.n = 1 << _p.level;
  for (; _p.level >= 0; _p.n /= 2, _p.level--)
    OMP_PARALLEL()
    Point point = _p;
    int _k;
    OMP(omp for schedule(static))
    for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
      point.i = _k;
      POINT_VARIABLES
@
@define end_foreach_fine_to_coarse() } OMP_END_PARALLEL() }

@def foreach_child() {
  int _i = 2*point.i - GHOSTS;
  point.level++;
  point.n *= 2;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2;
  point.level--;
  point.n /= 2;
}
@

#if TRASH
# undef trash
# define trash multigrid_trash
#endif

void multigrid_trash (void * alist)
{
  scalar * list = alist;
  Point p = *((Point *)grid);
  p.level = p.depth; p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < (p.n + 2*GHOSTS); i++)
      for (scalar s in list)
	((double *)(&p.d[p.level][i*datasize]))[s] = undefined;
}

void free_grid (void)
{
  if (!grid)
    return;
  Point * m = grid;
  for (int l = 0; l <= m->depth; l++)
    free (m->d[l]);
  free (m->d);
  free (m);
  grid = NULL;
}

void init_grid (int n)
{
  Point * m = grid;
  if (m && n == 1 << m->depth)
    return;
  free_grid();
  int r = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "multigrid: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    r++;
  }
  m = malloc(sizeof(Point));
  m->depth = r;
  N = 1 << r;
  m->d = malloc(sizeof(Point *)*(r + 1));
  for (int l = 0; l <= r; l++) {
    size_t len = _size(l)*datasize;
    m->d[l] = malloc (len);
    /* trash the data just to make sure it's either explicitly
       initialised or never touched */
    double * v = (double *) m->d[l];
    for (int i = 0; i < len/sizeof(double); i++)
      v[i] = undefined;
  }
  grid = m;
  trash (all);
  init_events();
}

void realloc_scalar (void)
{
  Point * p = grid;
  size_t oldatasize = datasize - sizeof(double);
  for (int l = 0; l <= p->depth; l++) {
    size_t len = _size(l);
    p->d[l] = realloc (p->d[l], len*datasize);
    char * data = p->d[l] + (len - 1)*oldatasize;
    for (int i = len - 1; i > 0; i--, data -= oldatasize)
      memmove (data + i*sizeof(double), data, oldatasize);  
  }
}

Point locate (double xp, double yp)
{
  Point point = *((Point *)grid);
  point.n = 1 << point.depth;
  double a = (xp - X0)/L0*point.n;
  point.i = a + GHOSTS;
  point.level = 
    (a >= 0.5 - GHOSTS && a < point.n + GHOSTS - 0.5) ? point.depth : - 1;
  return point;
}

#include "multigrid-common.h"

void multigrid1D_methods()
{
  multigrid_methods();
}
