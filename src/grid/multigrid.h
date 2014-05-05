#define GRIDNAME "Multigrid"

#define I      (point.i - GHOSTS)
#define J      (point.j - GHOSTS)
#define DELTA  (1./point.n)

struct _Point {
  char ** d;
  int level, depth;
  int i, j, n;
};
#define NN point.n // for output_stencil()

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
@def data(k,l)
  ((double *)&point.d[point.level][((point.i + k)*(point.n + 2*GHOSTS) +
				    point.j + l)*datasize]) @
@define allocated(...) true

/***** Multigrid variables and macros *****/
@define depth()       (((Point *)grid)->depth)
@def fine(a,k,l)
  ((double *)
   &point.d[point.level+1][((2*point.i-GHOSTS+k)*2*(point.n + GHOSTS) +
			    (2*point.j-GHOSTS+l))*datasize])[a]
@
@def coarse(a,k,l)
  ((double *)
   &point.d[point.level-1][(((point.i+GHOSTS)/2+k)*(point.n/2+2*GHOSTS) +
			    (point.j+GHOSTS)/2+l)*datasize])[a]
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1
  }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2; parent.j = (point.j + GHOSTS)/2;
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
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {
      POINT_VARIABLES
@
@define end_foreach_level() }} OMP_END_PARALLEL()

@def foreach(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = point.depth; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
    for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {
      POINT_VARIABLES
@
@define end_foreach() }} OMP_END_PARALLEL()

@def foreach_face_generic()
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  point.level = point.depth; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
    point.i = _k;
    for (point.j = GHOSTS; point.j <= point.n + GHOSTS; point.j++) {
      POINT_VARIABLES
@
@define end_foreach_face_generic() }} OMP_END_PARALLEL()

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

@define is_face_x() (point.j <= point.n)
@define is_face_y() (point.i <= point.n)

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
      for (point.j = GHOSTS; point.j < point.n + GHOSTS; point.j++) {
        POINT_VARIABLES
@
@define end_foreach_fine_to_coarse() }} OMP_END_PARALLEL() }

@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  point.n *= 2;
  for (int _k = 0; _k < 2; _k++)
    for (int _l = 0; _l < 2; _l++) {
      point.i = _i + _k; point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
  point.n /= 2;
}
@

@if TRASH
@ undef trash
@ define trash multigrid_trash
@endif

void multigrid_trash (void * alist)
{
  scalar * list = alist;
  Point p = *((Point *)grid);
  p.level = p.depth; p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < (p.n + 2*GHOSTS)*(p.n + 2*GHOSTS); i++)
      for (scalar s in list)
	((double *)(&p.d[p.level][i*datasize]))[s] = undefined;
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL()
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
  int _start = GHOSTS, _end = point.n + GHOSTS, _k;
  /* traverse corners only for top and bottom */
  if (d > left) { _start -= GHOSTS; _end += GHOSTS; }
  OMP(omp for schedule(static))
  for (_k = _start; _k < _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
    for (scalar s in list)
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
  }
  OMP_END_PARALLEL()
}

static void box_boundary_normal (const Boundary * b, vector * list)
{
  int d = ((BoxBoundary *)b)->d;
  int component = d/2; // index of normal component

  OMP_PARALLEL()
  Point point = *((Point *)grid);
  // set the indices of the normal face component in direction d
  if (d % 2) {
    ig = jg = 0;
  }
  else {
    ig = _ig[d]; jg = _jg[d];
  }
  point.level = depth(); point.n = 1 << point.level;
  int _start = GHOSTS, _end = point.n + GHOSTS, _k;
  OMP(omp for schedule(static))
  for (_k = _start; _k < _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
    for (vector v in list) {
      scalar s = (&v.x)[component];
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
    }
  }
  OMP_END_PARALLEL()
}

static void box_boundary_tangent (const Boundary * b, vector * list)
{
  int d = ((BoxBoundary *)b)->d;
  int component = (d/2 + 1) % 2; // index of tangential component

  OMP_PARALLEL()
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  point.level = depth(); point.n = 1 << point.level;
  int _start = GHOSTS, _end = point.n + GHOSTS, _k;
  OMP(omp for schedule(static))
  for (_k = _start; _k <= _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
    for (vector v in list) {
      scalar s = (&v.x)[component];
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
    }
  }
  OMP_END_PARALLEL()
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
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
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = b->restriction = box_boundary_level;
    b->normal  = box_boundary_normal;
    b->tangent = box_boundary_tangent;
    add_boundary (b);
  }  
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
  double b = (yp - Y0)/L0*point.n;
  point.j = b + GHOSTS;
  point.level = 
    (a >= 0.5 - GHOSTS && a < point.n + GHOSTS - 0.5 &&
     b >= 0.5 - GHOSTS && b < point.n + GHOSTS - 0.5) ? point.depth : - 1;
  return point;
}

#include "multigrid-common.h"
