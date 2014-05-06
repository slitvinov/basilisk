#define GRIDNAME "Cartesian 1D"
#define GHOSTS 1

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

@def foreach_face_generic(clause)
  OMP_PARALLEL()
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point = *((Point *)grid);
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 1; _k <= point.n + 1; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach_face_generic() } OMP_END_PARALLEL()

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

@define is_face_x() (true)
@define is_face_y() (point.i <= point.n)

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  point.i = d == right ? point.n : 1;
  {
    POINT_VARIABLES;
    for (scalar s in list)
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
  }
}

static void box_boundary_normal (const Boundary * b, vector * list)
{
  int d = ((BoxBoundary *)b)->d;
  int component = d/2; // index of normal component

  Point point = *((Point *)grid);
  // set the indices of the normal face component in direction d
  if (d % 2) {
    ig = jg = 0;
  }
  else {
    ig = _ig[d]; jg = _jg[d];
  }
  point.i = d == right ? point.n : 1;
  {
    POINT_VARIABLES;
    for (vector v in list) {
      scalar s = (&v.x)[component];
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
    }
  }
}

static void box_boundary_tangent (const Boundary * b, vector * list)
{
  // fixme: should not need this
  int d = ((BoxBoundary *)b)->d;
  int component = (d/2 + 1) % 2; // index of tangential component
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  point.i = d == right ? point.n : 1;
  {
    POINT_VARIABLES;
    for (vector v in list) {
      scalar s = (&v.x)[component];
      val(s,ig,jg) = _method[s].boundary[d] (point, s);
    }
  }
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
  if (p && p->n == n)
    return;
  free_grid();
  p = malloc(sizeof(Point));
  size_t len = (n + 2)*datasize;
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
  size_t len = (p->n + 2)*datasize;
  p->data = realloc (p->data, len);
  int oldatasize = datasize - sizeof(double);
  char * data = p->data + (p->n + 1)*oldatasize;
  for (int i = p->n + 1; i > 0; i--, data -= oldatasize)
    memmove (data + i*sizeof(double), data, oldatasize);
}

@if TRASH
@ undef trash
@ define trash cartesian1D_trash
@endif

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
  double a = (xp - X0)/L0*point.n;
  point.i = a + 1;
  point.level = (a > -0.5 && a < point.n + 0.5) ? 0 : - 1;
  return point;
}

void cartesian1D_methods()
{
  cartesian_methods();
}
