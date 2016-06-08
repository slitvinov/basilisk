#define GRIDNAME "Cartesian 1D"
#define dimension 1
#define GHOSTS 1

#define _I     (point.i - 1)
#define _DELTA (1./point.n)

typedef struct {
  Grid g;
  char * d;
  int n;
} Cartesian;

struct _Point {
  int i, j, level, n;
};
static Point last_point;

#define cartesian ((Cartesian *)grid)

@define data(k,l,m) ((double *)&cartesian->d[(point.i + k)*datasize])
@define allocated(...) true

@define POINT_VARIABLES VARIABLES

@def foreach()
  OMP_PARALLEL() {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 1; _k <= point.n; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach() }}

@def foreach_face_generic()
  OMP_PARALLEL() {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 1; _k <= point.n + 1; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach_face_generic() }}

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

@define is_face_x() (true)
@define is_face_y() (point.i <= point.n)

// ghost cell coordinates for each direction
static int _ig[] = {1,-1};

// Box boundaries

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  if (!list)
    return;
  
  int d = ((BoxBoundary *)b)->d;

  Point point;
  point.n = cartesian->n;
  ig = _ig[d];
  assert (d <= left);
  point.i = d == right ? point.n + GHOSTS : GHOSTS;
  Point neighbor = {point.i + ig};
  for (scalar s in list) {
    scalar b = s.v.x;
    val(s,ig) = b.boundary[d] (point, neighbor, s);
  }
}

static double periodic_bc (Point point, Point neighbor, scalar s);

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[d] != periodic_bc) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  if (centered) {
    Point point;
    point.n = cartesian->n;
    ig = _ig[d];
    point.i = d == right ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig};
    for (scalar s in centered)
      val(s,ig) = s.boundary[d] (point, neighbor, s);
    free (centered);
  }
    
  box_boundary_level_normal (b, normal, l);
  free (normal);
}

// periodic boundaries

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  Point point = *((Point *)grid);
  point.i = 0;
  for (int i = 0; i < GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i + point.n];
  for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i - point.n];

  free (list1);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  free (cartesian->d);
  free (cartesian);
  grid = NULL;
}

void init_grid (int n)
{
  if (cartesian && n == cartesian->n)
    return;
  free_grid();
  Cartesian * p = malloc(sizeof(Cartesian));
  size_t len = (n + 2)*datasize;
  p->n = N = n;
  p->d = malloc (len);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->d;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = (Grid *) p;
  trash (all);
  // box boundaries
  for (int d = 0; d < 2; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
    add_boundary (b);
  }
  // periodic boundaries
  Boundary * b = calloc (1, sizeof (Boundary));
  b->level = periodic_boundary_level_x;
  add_boundary (b);
  // mesh size
  grid->n = grid->tn = n;
}

void realloc_scalar (void)
{
  Cartesian * p = cartesian;
  size_t len = (p->n + 2)*datasize;
  p->d = realloc (p->d, len);
  int oldatasize = datasize - sizeof(double);
  char * data = p->d + (p->n + 1)*oldatasize;
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
  char * data = cartesian->d;
  for (int i = 0; i < cartesian->n + 2; i++, data += datasize) {
    double * v = (double *) data;
    for (scalar s in list)
      if (!is_constant(s))
	v[s.i] = undefined;
  }
}

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point;
  point.n = cartesian->n;
  double a = (p.x - X0)/L0*point.n;
  point.i = a + 1;
  point.level = (a > -0.5 && a < point.n + 0.5) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"

void cartesian1D_methods()
{
  cartesian_methods();
}
