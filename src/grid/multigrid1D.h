#define GRIDNAME "Multigrid 1D"
#define dimension 1
#define GHOSTS 1

#define I      (point.i - GHOSTS)
#define DELTA  (1./point.n)

typedef struct {
  char ** d;
  int depth;
} Multigrid;

struct _Point {
  int i, level, n;
};
static Point last_point;

#define multigrid ((Multigrid *)grid)

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*datasize]))

/***** Cartesian macros *****/
@def data(k,l,m)
((double *)&multigrid->d[point.level][(point.i + k)*datasize + (l) - (l)]) @
@define allocated(...) true
@define allocated_child(...) true

/***** Multigrid variables and macros *****/
@define depth()       (((Multigrid *)grid)->depth)
@def fine(a,k,l,m)
  ((double *)
   &multigrid->d[point.level+1][(2*point.i-GHOSTS+k)*datasize])[a.i]
@
  @def coarse(a,k,l,m)
  ((double *)
   &multigrid->d[point.level-1][((point.i+GHOSTS)/2+k)*datasize])[a.i]
@
@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x; } child = { 2*((point.i+GHOSTS)%2)-1 }; NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
@
@def foreach_level(l)
  OMP_PARALLEL()
  int ig = 0; NOT_UNUSED(ig);
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
  int ig = 0; NOT_UNUSED(ig);
  Point point = *((Point *)grid);
  point.level = multigrid->depth; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach() } OMP_END_PARALLEL()

@def foreach_face_generic(clause)
  OMP_PARALLEL()
  int ig = 0; NOT_UNUSED(ig);
  Point point = *((Point *)grid);
  point.level = multigrid->depth; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
    point.i = _k;
    POINT_VARIABLES
@
@define end_foreach_face_generic() } OMP_END_PARALLEL()

@define is_face_x() true

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

@define is_coarse() (point.level < depth())

@def foreach_fine_to_coarse() {
  int ig = 0; NOT_UNUSED(ig);
  Point _p;
  _p.level = multigrid->depth - 1; _p.n = 1 << _p.level;
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
  Point p;
  p.level = multigrid->depth; p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < (p.n + 2*GHOSTS); i++)
      for (scalar s in list)
	if (!is_constant(s))
	  ((double *)(&multigrid->d[p.level][i*datasize]))[s.i] = undefined;
}

// ghost cell coordinates for each direction
static int _ig[] = {1,-1};

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  Point point = *((Point *)grid);
  ig = d % 2 ? 0 : _ig[d];
  point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
  assert (d <= left);
  point.i = d == right ? point.n + GHOSTS : GHOSTS;
  Point neighbor = {point.i + ig, point.level};
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
    Point point = *((Point *)grid);
    ig = _ig[d];
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
    assert (d <= left);
    point.i = d == right ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig, point.level};
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

  Point point = {0};
  point.level = l < 0 ? depth() : l; point.n = 1 << point.level;
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
  Multigrid * m = grid;
  for (int l = 0; l <= m->depth; l++)
    free (m->d[l]);
  free (m->d);
  free (m);
  grid = NULL;
}

void init_grid (int n)
{
  Multigrid * m = grid;
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
  m = malloc(sizeof(Multigrid));
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
  // box boundaries
  for (int d = 0; d <= left; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level = b->restriction = box_boundary_level;
    add_boundary (b);
  }
  // periodic boundaries
  Boundary * b = calloc (1, sizeof (Boundary));
  b->level = b->restriction = periodic_boundary_level_x;
  add_boundary (b);
  init_events();
}

void realloc_scalar (void)
{
  Multigrid * p = grid;
  size_t oldatasize = datasize - sizeof(double);
  for (int l = 0; l <= p->depth; l++) {
    size_t len = _size(l);
    p->d[l] = realloc (p->d[l], len*datasize);
    char * data = p->d[l] + (len - 1)*oldatasize;
    for (int i = len - 1; i > 0; i--, data -= oldatasize)
      memmove (data + i*sizeof(double), data, oldatasize);  
  }
}

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point = *((Point *)grid);
  point.n = 1 << multigrid->depth;
  double a = (p.x - X0)/L0*point.n;
  point.i = a + GHOSTS;
  point.level = 
    (a >= 0.5 - GHOSTS && a < point.n + GHOSTS - 0.5) ? multigrid->depth : - 1;
  return point;
}

#include "multigrid-common.h"

void multigrid1D_methods()
{
  multigrid_methods();
}
