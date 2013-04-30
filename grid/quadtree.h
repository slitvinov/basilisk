#define GRIDNAME "Quadtree"

#define TWO_ONE 1 // enforce 2:1 refinement ratio

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

typedef struct {
  char flags, neighbors; // number of neighboring leaves
} Cell;

enum {
  active = 1 << 0,
  leaf   = 1 << 1
};

#define _CORNER 4
#define is_leaf(cell)    ((cell).flags & leaf)
#define is_active(cell)  ((cell).flags & active)
#define is_refined(cell) (is_active(cell) && !is_leaf(cell))
#define is_corner(cell)  (stage == _CORNER)

typedef struct _Quadtree Point;
typedef struct _Quadtree Quadtree;
typedef struct _Cache    Cache;

struct _Quadtree {
  int depth;        /* the maximum depth of the tree */

  Quadtree * back;  /* back pointer to the "parent" quadtree */
  char ** m;        /* the grids at each level */
  int i, j, level;  /* the current cell index and level */

  Cache * halo;     /* halo indices for each level */
  Cache * active;   /* active cells indices for each level */

  bool dirty;       /* whether caches should be updated */
};

typedef struct {
  int i, j;
} Index;

struct _Cache {
  Index * p;
  int n, nm;
};

static void cache_append (Cache * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (Index)*c->nm);
  }
  c->p[c->n].i = p.i;
  c->p[c->n].j = p.j;
  c->n++;
}

static void cache_init (Cache * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))

/***** Multigrid macros *****/
#define depth()      (((Quadtree *)grid)->depth)
#define aparent(k,l) \
  CELL(point.m, point.level-1, ((point.i+GHOSTS)/2+k)*(_n/2+2*GHOSTS) + \
       (point.j+GHOSTS)/2+l)
#define child(k,l)   \
  CELL(point.m, point.level+1, (2*point.i-GHOSTS+k)*2*(_n + GHOSTS) +	\
       (2*point.j-GHOSTS+l))

/***** Quadtree macros ****/
#define _n (1 << point.level) /* fixme */
#define cell							\
  CELL(point.m, point.level, point.i*(_n + 2*GHOSTS) + point.j)
#define _neighbor(k,l)							\
  CELL(point.m, point.level, (point.i + k)*(_n + 2*GHOSTS) + point.j + l)
#define parent             aparent(0,0)

/***** Data macros *****/
#define data(k,l)							\
  ((double *) &point.m[point.level][((point.i + k)*(_n + 2*GHOSTS) +	\
				     (point.j + l))*(sizeof(Cell) + datasize) \
				    + sizeof(Cell)])
#define field(cell) ((double *)(((char *) &cell) + sizeof(Cell)))
#define _fine(a,k,l) field(child(k,l))[a]
#define _coarse(a,k,l) field(aparent(k,l))[a]

#define POINT_VARIABLES						     \
  VARIABLES							     \
  int level = point.level; NOT_UNUSED(level);			     \
  struct { int x, y; } child = {				     \
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1		     \
  }; NOT_UNUSED(child);

/* ===============================================================
 *                    Quadtree traversal
 * recursive() below is for reference only. The macro
 * foreach_cell() is a stack-based implementation of
 * recursive(). It is about 12% slower than the recursive
 * version and 60% slower than simple array traversal.
 *
 * This article was useful:
 * http://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
 *
 * =============================================================== */

#define _BOTTOM (2*point.j - GHOSTS)
#define _TOP    (_BOTTOM + 1)
#define _LEFT   (2*point.i - GHOSTS)
#define _RIGHT  (_LEFT + 1)

void recursive (Point point)
{
  if (point.level == point.depth) {
    /* do something */
  }
  else {
    Point p1 = point; p1.level = point.level + 1;
    p1.i = _LEFT;  p1.j = _TOP;    recursive (p1);
    p1.i = _RIGHT; p1.j = _TOP;    recursive (p1);
    p1.i = _LEFT;  p1.j = _BOTTOM; recursive (p1);
    p1.i = _RIGHT; p1.j = _BOTTOM; recursive (p1);
  }
}

#define STACKSIZE 20
#define _push(b,c,d,e)					                \
  { _s++; stack[_s].l = b; stack[_s].i = c; stack[_s].j = d;		\
    stack[_s].stage = e; }
#define _pop(b,c,d,e)							\
  { b = stack[_s].l; c = stack[_s].i; d = stack[_s].j;			\
    e = stack[_s].stage; _s--; }

#define foreach_cell()							\
  {									\
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);			\
    Quadtree point = *((Quadtree *)grid); point.back = grid;		\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;	\
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
        POINT_VARIABLES;						\
	/* do something */
#define end_foreach_cell()						\
        if (point.level < point.depth) {				\
	  _push (point.level, point.i, point.j, 1);			\
          _push (point.level + 1, _LEFT, _TOP, 0);			\
        }								\
        break;								\
      }									\
      case 1: _push (point.level, point.i, point.j, 2);			\
              _push (point.level + 1, _RIGHT, _TOP,    0); break;	\
      case 2: _push (point.level, point.i, point.j, 3);		        \
              _push (point.level + 1, _LEFT,  _BOTTOM, 0); break;	\
      case 3: _push (point.level + 1, _RIGHT, _BOTTOM, 0); break;	\
      }								        \
    }                                                                   \
  }

#define foreach_cell_post(condition)					\
  {									\
    Quadtree point = *((Quadtree *)grid); point.back = grid;		\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;	\
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
        POINT_VARIABLES;						\
	if (condition) {						\
	  if (point.level == point.depth)	{			\
	    _push (point.level, point.i, point.j, 4);			\
	  }								\
	  else {							\
	    _push (point.level, point.i, point.j, 1);			\
	    _push (point.level + 1, _LEFT, _TOP, 0);			\
	  }								\
	}								\
	break;								\
      }									\
      case 1: _push (point.level, point.i, point.j, 2);                 \
              _push (point.level + 1, _RIGHT, _TOP,    0); break;	\
      case 2: _push (point.level, point.i, point.j, 3);                 \
	      _push (point.level + 1, _LEFT,  _BOTTOM, 0); break;	\
      case 3: _push (point.level, point.i, point.j, 4);                 \
	      _push (point.level + 1, _RIGHT, _BOTTOM, 0); break;	\
      case 4: {								\
        POINT_VARIABLES;						\
	/* do something */
#define end_foreach_cell_post()						\
      }									\
      }								        \
    }                                                                   \
  }

#define foreach_boundary_cell(dir)					\
  {									\
    int ig = _ig[dir], jg = _jg[dir];	NOT_UNUSED(ig); NOT_UNUSED(jg);	\
    Quadtree point = *((Quadtree *)grid); point.back = grid;		\
    int _d = dir; NOT_UNUSED(_d);					\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;	\
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0:								\
        /* corners */							\
        if (_d < top) {							\
  	  if (point.j == GHOSTS)					\
	    _push (point.level, point.i, point.j - 1, _CORNER);		\
	  if (point.j == _n + 2*GHOSTS - 2)			        \
	    _push (point.level, point.i, point.j + 1, _CORNER);		\
	} else {							\
	  if (point.i == GHOSTS)					\
	    _push (point.level, point.i - 1, point.j, _CORNER);		\
	  if (point.i == _n + 2*GHOSTS - 2)			        \
	    _push (point.level, point.i + 1, point.j, _CORNER);		\
        }							        \
	/* fall through */						\
      case _CORNER: {							\
          POINT_VARIABLES;						\
  	  /* do something */
#define end_foreach_boundary_cell()					\
        }								\
	if (stage == _CORNER) continue;				        \
        /* children */							\
        if (point.level < point.depth) {                                \
	  _push (point.level, point.i, point.j, 1);			\
	  int k = _d > left ? _LEFT : _RIGHT - _d;			\
	  int l = _d < top  ? _TOP  : _TOP + 2 - _d;			\
	  _push (point.level + 1, k, l, 0);				\
	}								\
	break;								\
      case 1: {								\
  	  int k = _d > left ? _RIGHT : _RIGHT - _d;			\
	  int l = _d < top  ? _BOTTOM  : _TOP + 2 - _d;			\
	  _push (point.level + 1, k, l, 0);				\
	  break;							\
        }								\
      }									\
    }                                                                   \
  }

#define foreach(clause)     {						\
  update_cache();							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  OMP_PARALLEL()							\
  Quadtree point = *((Quadtree *)grid); point.back = grid;		\
  for (int _l = 0; _l <= depth(); _l++)					\
    OMP(omp for schedule(static) clause)				\
    for (int _k = 0; _k < point.active[_l].n; _k++) {			\
      point.i = point.active[_l].p[_k].i;				\
      point.j = point.active[_l].p[_k].j;				\
      point.level = _l;							\
      POINT_VARIABLES;							\
      if (is_leaf (cell)) {
#define end_foreach() } } OMP_END_PARALLEL() }

#define foreach_fine_to_coarse(clause)     {				\
  update_cache();							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  OMP_PARALLEL()							\
  Quadtree point = *((Quadtree *)grid); point.back = grid;		\
  for (int _l = depth() - 1; _l >= 0; _l--)				\
    OMP(omp for schedule(static) clause)				\
    for (int _k = 0; _k < point.active[_l].n; _k++) {			\
      point.i = point.active[_l].p[_k].i;			        \
      point.j = point.active[_l].p[_k].j;				\
      point.level = _l;							\
      POINT_VARIABLES;							\
      if (!is_leaf (cell)) {
#define end_foreach_fine_to_coarse() } } OMP_END_PARALLEL() }

#define foreach_level(l)     {						\
  update_cache();							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  OMP_PARALLEL()							\
  Quadtree point = *((Quadtree *)grid); point.back = grid;		\
  for (int _l = l; _l >= 0; _l--)					\
    OMP(omp for schedule(static))					\
    for (int _k = 0; _k < point.active[_l].n; _k++) {			\
      point.i = point.active[_l].p[_k].i;				\
      point.j = point.active[_l].p[_k].j;				\
      point.level = _l;							\
      POINT_VARIABLES;							\
      if (_l == l || is_leaf (cell)) {
#define end_foreach_level() } } OMP_END_PARALLEL() }

#define foreach_leaf()            foreach_cell() if (is_leaf (cell)) {
#define end_foreach_leaf()        continue; } end_foreach_cell()

void alloc_layer (Quadtree * p)
{
  Quadtree * q = p->back;
  q->depth++; p->depth++;
  q->m = &(q->m[-1]);
  q->m = realloc(q->m, sizeof (char *)*(q->depth + 2)); 
  q->m = &(q->m[1]);
  p->m = q->m;
  q->m[q->depth] = calloc (_size(q->depth), sizeof (Cell) + datasize);
  q->active = realloc (q->active, (q->depth + 1)*sizeof (Cache));
  cache_init (&q->active[q->depth]);
  q->halo = realloc (q->halo, (q->depth + 1)*sizeof (Cache));
  cache_init (&q->halo[q->depth]);
#if TRASH
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  foreach_cell()
    if (level == q->depth) {
      for (scalar v = 0; v < nvar; v++)
	val(v,0,0) = undefined;
      continue;
    }
#endif
}

Point refine_cell (Point point, scalar * list)
{
#if TWO_ONE
  /* refine neighborhood if required */
  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)
      for (int l = 0; l != 2*child.y; l += child.y)
	if (aparent(k,l).flags & leaf) {
	  Point p = point;
	  /* fixme: this should be made
	     independent from the quadtree implementation */
	  p.level = point.level - 1;
	  p.i = (point.i + GHOSTS)/2 + k;
	  p.j = (point.j + GHOSTS)/2 + l;
	  p = refine_cell (p, list);
	  assert (p.m == point.m);
	}
#endif

  /* refine */
  point.back->dirty = true;
  if (point.level == point.depth) alloc_layer(&point);
  cell.flags &= ~leaf;
  /* update neighborhood */
  for (int o = -GHOSTS; o <= GHOSTS; o++)
    for (int p = -GHOSTS; p <= GHOSTS; p++)
      neighbor(o,p).neighbors--;

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      assert(!(child(k,l).flags & active));
      child(k,l).flags |= (active | leaf);
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors++;
#if 0
      /* bilinear interpolation from coarser level */
      for (scalar v in list)
	fine(v,k,l) = 
	  (9.*v[] + 3.*(v[2*k-1,0] + v[0,2*l-1]) + v[2*k-1,2*l-1])/16.;
#else
      /* linear interpolation from coarser level (conservative) */
      for (scalar v in list)
	fine(v,k,l) = v[] + ((v[1,0] - v[-1,0])*(2*k-1)/8. +
			     (v[0,1] - v[0,-1])*(2*l-1)/8.);
#endif
    }

  return point;
}

bool coarsen_cell (Point point)
{
#if TWO_ONE
  /* check that neighboring cells are not too fine */
  for (int k = -1; k < 3; k++)
    for (int l = -1; l < 3; l++)
      if (is_active (child(k,l)) && !is_leaf (child(k,l)))
	return false; // cannot coarsen
#endif

  /* coarsen */
  point.back->dirty = true;
  cell.flags |= leaf;
  /* update neighborhood */
  for (int o = -GHOSTS; o <= GHOSTS; o++)
    for (int p = -GHOSTS; p <= GHOSTS; p++)
      neighbor(o,p).neighbors++;

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      child(k,l).flags &= ~(leaf|active);
#if TRASH
      /* trash the data just to make sure it's never touched */
      for (scalar v = 0; v < nvar; v++)
	fine(v,k,l) = undefined;
#endif
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors--;
    }
  return true;
}

static void update_cache (void)
{
  Quadtree * q = grid;
  if (!q->dirty)
    return;

  /* empty caches */
  for (int l = 0; l <= depth(); l++)
    q->halo[l].n = q->active[l].n = 0;

  foreach_cell() {
    if (!is_active (cell)) {
      if (cell.neighbors > 0)
	/* update halo cache (prolongation) */
	cache_append (&q->halo[level], point);
      else
	continue;
    }
    else {
      if (!is_leaf (cell) && cell.neighbors > 0)
	/* update halo cache (restriction) */
	cache_append (&q->halo[level], point);
      /* update active cache */
      cache_append (&q->active[level], point);
    }
  }

  q->dirty = false;
}

/* breadth-first traversal of halos from coarse to fine */
#define foreach_halo_coarse_to_fine(depth1)    {			\
  update_cache();							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  int _depth = depth1 < 0 ? depth() : depth1;				\
  OMP_PARALLEL()							\
  Quadtree point = *((Quadtree *)grid); point.back = grid;		\
  for (int _l = 0; _l <= _depth; _l++)                                  \
    OMP(omp for schedule(static))					\
    for (int _k = 0; _k < point.halo[_l].n; _k++) {			\
      point.i = point.halo[_l].p[_k].i;					\
      point.j = point.halo[_l].p[_k].j;					\
      point.level = _l;							\
      POINT_VARIABLES;							\
      if (!is_active(cell)) {
#define end_foreach_halo_coarse_to_fine()	\
  } } OMP_END_PARALLEL() }

/* breadth-first traversal of halos from fine to coarse */
#define foreach_halo_fine_to_coarse()    {				\
  update_cache();							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  OMP_PARALLEL()							\
  Quadtree point = *((Quadtree *)grid); point.back = grid;		\
  for (int _l = depth() - 1; _l >= 0; _l--)				\
    OMP(omp for schedule(static))					\
    for (int _k = 0; _k < point.halo[_l].n; _k++) {			\
      point.i = point.halo[_l].p[_k].i;					\
      point.j = point.halo[_l].p[_k].j;					\
      point.level = _l;							\
      POINT_VARIABLES;							\
      if (is_active(cell)) {
#define end_foreach_halo_fine_to_coarse()	\
  } } OMP_END_PARALLEL() }

void init_grid (int n)
{
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Quadtree * q = malloc (sizeof (Quadtree));
  q->depth = 0; q->i = q->j = GHOSTS; q->level = 0.;
  q->m = malloc(sizeof (char *)*2);
  /* make sure we don't try to access level -1 */
  q->m[0] = NULL; q->m = &(q->m[1]);
  /* initialise the root cell */
  q->m[0] = calloc (_size(0), sizeof (Cell) + datasize);
  CELL(q->m, 0, 2 + 2*GHOSTS).flags |= (leaf | active);
  CELL(q->m, 0, 2 + 2*GHOSTS).neighbors = 1; // only itself as neighbor
  q->active = calloc (1, sizeof (Cache));
  q->halo = calloc (1, sizeof (Cache));
  q->dirty = true;
  grid = q;
  while (depth--)
    foreach_leaf()
      point = refine_cell (point, all);
  update_cache();
  init_solver();
}

void free_grid (void)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++) {
    free (q->m[l]);
    free (q->halo[l].p);
    free (q->active[l].p);
  }
  q->m = &(q->m[-1]);
  free(q->m);
  free(q->halo);
  free(q->active);
  free(q);
  free_boundaries();
}

void output_cells (FILE * fp);

void check_two_one (void)
{
  foreach_leaf()
    if (level > 0)
      for (int k = -1; k <= 1; k++)
	for (int l = -1; l <= 1; l++) {
	  /* fixme: all this mess is just to ignore ghost cells */
	  int i = (point.i + GHOSTS)/2 + k;
	  int j = (point.j + GHOSTS)/2 + l;
	  double x = ((i - GHOSTS + 0.5)*DELTA*2. - 0.5);
	  double y = ((j - GHOSTS + 0.5)*DELTA*2. - 0.5);
	  if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 && 
	      !(aparent(k,l).flags & active)) {
	    FILE * fp = fopen("check_two_one_loc", "w");
	    fprintf (fp,
		     "# %d %d\n"
		     "%g %g\n%g %g\n",
		     k, l,
		     ((I + 0.5)*DELTA - 0.5),
		     ((J + 0.5)*DELTA - 0.5),
		     x, y);
	    fclose (fp);
#if 0
	    fp = fopen("check_two_one", "w");
	    output_cells (fp);
	    fclose (fp);
#endif
	    assert (false);
	  }
	}
}

#include "quadtree-common.h"
