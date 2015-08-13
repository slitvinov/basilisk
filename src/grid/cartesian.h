#define GRIDNAME "Cartesian"
#define dimension 2
#define GHOSTS 1

#define I     (point.i - 1)
#define J     (point.j - 1)
#define DELTA (1./point.n)

struct _Point {
  int i, j;
  int level; // only to return level in locate()
  int n;
  char * data;
};
static Point last_point;

@def data(k,l,m) ((double *)&point.data[((point.i + k)*(point.n + 2) +
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
  for (int i = 0; i < sq(p->n + 2); i++)
    for (scalar s in list)
      if (!is_constant(s))
	((double *)(&p->data[i*datasize]))[s.i] = undefined;
}

// ghost cell coordinates for each direction
static int _ig[] = {1,-1,0,0}, _jg[] = {0,0,1,-1};

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL();
  Point point = *((Point *)grid);
  if (d % 2)
    ig = jg = 0;
  else {
    ig = _ig[d]; jg = _jg[d];
  }
  int _start = GHOSTS, _end = point.n + GHOSTS, _k;  
  OMP(omp for schedule(static))
  for (_k = _start; _k < _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig, point.j + jg};
    for (scalar s in list) {
      scalar b = s.v.x;
      val(s,ig,jg) = b.boundary[d] (point, neighbor, s);
    }
  }
  OMP_END_PARALLEL();
}

static void box_boundary_level_tangent (const Boundary * b, 
					scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL();
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  int _start = GHOSTS, _end = point.n + 2*GHOSTS, _k;
  
  OMP(omp for schedule(static))
  for (_k = _start; _k < _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
    point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig, point.j + jg};
    for (scalar s in list) {
      scalar b = s.v.y;
      val(s,ig,jg) = b.boundary[d] (point, neighbor, s);
    }
  }
  OMP_END_PARALLEL();
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL, * tangent = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
	else {
	  scalar b = s.v.y;
	  if (b.boundary[d])
	    tangent = list_add (tangent, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  OMP_PARALLEL();
  Point point = *((Point *)grid);
  ig = _ig[d]; jg = _jg[d];
  int _start = 1, _end = point.n, _k;
  /* traverse corners only for top and bottom */
  if (d > left) { _start--; _end++; }
  OMP(omp for schedule(static))
  for (_k = _start; _k <= _end; _k++) {
    point.i = d > left ? _k : d == right ? point.n : 1;
    point.j = d < top  ? _k : d == top   ? point.n : 1;
    Point neighbor = {point.i + ig, point.j + jg};
    for (scalar s in centered) {
      scalar b = (s.v.x.i < 0 ? s :
		  s.i == s.v.x.i && d < top ? s.v.x :
		  s.i == s.v.y.i && d >= top ? s.v.x :
		  s.v.y);
      val(s,ig,jg) = b.boundary[d] (point, neighbor, s);
    }
  }
  OMP_END_PARALLEL();
  free (centered);

  box_boundary_level_normal (b, normal, l);
  free (normal);
  box_boundary_level_tangent (b, tangent, l);
  free (tangent);
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
    add_boundary (b);
  }
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

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point = *((Point *)grid);
  point.i = (p.x - X0)/L0*point.n + 1;
  point.j = (p.y - Y0)/L0*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n &&
		 point.j >= 1 && point.j <= point.n) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"
