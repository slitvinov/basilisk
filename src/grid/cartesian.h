#define GRIDNAME "Cartesian"

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

struct _Point {
  char * data;
  int i, j, n;
  int level; // only to return level in locate()
};
#define NN point.n // for output_stencil()

@def data(k,l) ((double *)&point.data[((point.i + k)*(point.n + 2) +
				       (point.j + l))*datasize]) @
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
    for (point.j = 1; point.j <= point.n; point.j++) {
      POINT_VARIABLES
@
@define end_foreach() }} OMP_END_PARALLEL()

@def foreach_face_generic(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 1; _k <= point.n + 1; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n + 1; point.j++) {
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

@if TRASH
@ undef trash
@ define trash cartesian_trash
@endif

void cartesian_trash (void * alist)
{
  scalar * list = alist;
  Point * p = grid;
  for (int i = 0; i < (p->n + 2)*(p->n + 2); i++)
    for (scalar s in list)
      ((double *)(&p->data[i*datasize]))[s] = undefined;
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL()
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  int _start = 1, _end = point.n, _k;
  /* traverse corners only for top and bottom */
  if (d > left) { _start--; _end++; }
  OMP(omp for schedule(static))
  for (_k = _start; _k <= _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n : 1;
    point.j = d < top  ? _k : d == top   ? point.n : 1;
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
  int _start = 1, _end = point.n + 1, _k;
  OMP(omp for schedule(static))
  for (_k = _start; _k < _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n : 1;
    point.j = d < top  ? _k : d == top   ? point.n : 1;
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
  int _start = 1, _end = point.n + 1, _k;
  OMP(omp for schedule(static))
  for (_k = _start; _k <= _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n : 1;
    point.j = d < top  ? _k : d == top   ? point.n : 1;
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
  Point * p = grid;
  free (p->data);
  free (p);
  grid = NULL;
}

void init_grid (int n)
{
  Point * p = grid;
  if (p && n == p->n)
    return;
  free_grid();
  p = malloc(sizeof(Point));
  size_t len = (n + 2)*(n + 2)*datasize;
  p->n = N = n;
  p->data = malloc (len);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->data;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = p;
  trash (all);
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
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
  size_t len = (p->n + 2)*(p->n + 2);
  p->data = realloc (p->data, len*datasize);
  char * data = p->data + (len - 1)*oldatasize;
  for (int i = len - 1; i > 0; i--, data -= oldatasize)
    memmove (data + i*sizeof(double), data, oldatasize);  
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
