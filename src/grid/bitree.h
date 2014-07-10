#define GRIDNAME "Binary tree"
#define dimension 1

#define DYNAMIC 1 // use dynamic data allocation
#define TWO_ONE 1 // enforce 2:1 refinement ratio

#define I     (point.i - GHOSTS)
#define J     -0.5
#define DELTA (1./(1 << point.level))

typedef struct {
  int flags, neighbors; // number of neighboring leaves
} Cell;

enum {
  active  = 1 << 0,
  leaf    = 1 << 1,
  fghost  = 1 << 2,
  refined = 1 << 3
};

#define is_active(cell)  ((cell).flags & active)
#define is_leaf(cell)    ((cell).flags & leaf)
#define is_ghost(cell)   ((cell).flags & fghost)
#define is_refined(cell) (is_active(cell) && !is_leaf(cell))
// #define is_corner(cell)  (stage == _CORNER)
#define is_coarse()      (!is_leaf(cell))

// Caches

typedef struct {
  int i;
} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

typedef struct {
  int i, level;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;

static void cache_init (Cache * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

// Bitree

typedef struct _Point Bitree;
struct _Point {
  int depth;        /* the maximum depth of the tree */

  Bitree * back;  /* back pointer to the "parent" bitree */
#if DYNAMIC
  char *** m;       /* the grids at each level */
#else
  char ** m;        /* the grids at each level */
#endif
  int i, j, level;  /* the current cell index and level */

  Cache        leaves;  /* leaf indices */
  CacheLevel * halo;    /* halo indices for each level */
  CacheLevel * active;  /* active cells indices for each level */

  bool dirty;       /* whether caches should be updated */
};

static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (IndexLevel)*c->nm);
  }
  c->p[c->n].i = p.i;
  c->n++;
}

static void cache_append (Cache * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (Index)*c->nm);
  }
  c->p[c->n].i = p.i;
  c->p[c->n].level = p.level;
  c->n++;
}

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n;
}

#if DYNAMIC
  @define CELL(m,level,i) (*((Cell *)m[level][i]))
  @define allocated(i,j) (point.m[point.level][_index(i,j)])
#else
  @define CELL(m,level,i) (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))
  @define allocated(i,j) true
#endif

/***** Multigrid macros *****/
@define depth()      (((Bitree *)grid)->depth)
@define _index(k,l)  (point.i + k + (l) - (l))
@define _parentindex(k,l) ((point.i + GHOSTS)/2 + k + (l) - (l))
@define _childindex(k,l) (2*point.i - GHOSTS + k + (l) - (l))
@define aparent(k,l) CELL(point.m, point.level-1, _parentindex(k,l))
@define child(k,l)   CELL(point.m, point.level+1, _childindex(k,l))

/***** Bitree macros ****/
@define NN (1 << point.level)
@define cell		CELL(point.m, point.level, _index(0,0))
@define neighbor(k,l)	CELL(point.m, point.level, _index(k,l))

/***** Data macros *****/
#if DYNAMIC
  @def data(k,l)
    ((double *) (point.m[point.level][_index(k,l)] + sizeof(Cell))) @
  @def fine(a,k,l)
    ((double *) (point.m[point.level+1][_childindex(k,l)] + sizeof(Cell)))[a] @
  @def coarse(a,k,l)
    ((double *) (point.m[point.level-1][_parentindex(k,l)] + sizeof(Cell)))[a] @
#else // !DYNAMIC
  @def data(k,l)
    ((double *) &point.m[point.level]
     [_index(k,l)*(sizeof(Cell) + datasize) + sizeof(Cell)]) @
  @def fine(a,k,l)
    ((double *) &point.m[point.level+1]
     [_childindex(k,l)*(sizeof(Cell) + datasize) + sizeof(Cell)])[a] @
  @def coarse(a,k,l)
    ((double *) &point.m[point.level-1]
     [_parentindex(k,l)*(sizeof(Cell) + datasize) + sizeof(Cell)])[a] @
#endif // !DYNAMIC

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

/* ===============================================================
 *                    Bitree traversal
 * recursive() below is for reference only. The macro
 * foreach_cell() is a stack-based implementation of
 * recursive(). It is about 12% slower than the recursive
 * version and 60% slower than simple array traversal.
 *
 * This article was useful:
 * http://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
 *
 * =============================================================== */

#define _LEFT   (2*point.i - GHOSTS)
#define _RIGHT  (_LEFT + 1)

void recursive (Point point)
{
  if (point.level == point.depth) {
    /* do something */
  }
  else {
    Point p1 = point; p1.level = point.level + 1;
    p1.i = _LEFT;  recursive (p1);
    p1.i = _RIGHT; recursive (p1);
  }
}

#define STACKSIZE 20
#define _push(b,c,e)					\
  { _s++; stack[_s].l = b; stack[_s].i = c;		\
    stack[_s].stage = e; }
#define _pop(b,c,e)							\
  { b = stack[_s].l; c = stack[_s].i;					\
    e = stack[_s].stage; _s--; }

@def foreach_cell()
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Bitree point = *((Bitree *)grid); point.back = grid;
    struct { int l, i, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell()
        if (point.level < point.depth) {
	  _push (point.level, point.i, 1);
          _push (point.level + 1, _LEFT, 0);
        }
        break;
      }
      case 1: _push (point.level + 1, _RIGHT, 0); break;
      }
    }
  }
@

@def foreach_cell_post(condition)
  {
    Bitree point = *((Bitree *)grid); point.back = grid;
    struct { int l, i, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
	if (condition) {
	  if (point.level == point.depth) {
	    _push (point.level, point.i, 4);
	  }
	  else {
	    _push (point.level, point.i, 1);
	    _push (point.level + 1, _LEFT, 0);
	  }
	}
	break;
      }
      case 1: _push (point.level + 1, _RIGHT, 0); break;
      case 4: {
        POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell_post()
      }
      }
    }
  }
@

@define corners()
@def foreach_boundary_cell(dir,corners)
  if (dir <= left) { _OMPSTART /* for face reduction */
    int ig = _ig[dir], jg = _jg[dir];	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Bitree point = *((Bitree *)grid); point.back = grid;
    int _d = dir;
    struct { int l, i, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
          POINT_VARIABLES;
  	  /* do something */
@
@def end_foreach_boundary_cell()
        }
        /* children */
        if (point.level < point.depth)
	  _push (point.level + 1, _d == left ? _LEFT : _RIGHT, 0);
	break;
      }
    }  _OMPEND
  }
@

@def foreach_boundary_cell_post(dir,condition)
  if (dir <= left) {
    int ig = _ig[dir], jg = _jg[dir];	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Bitree point = *((Bitree *)grid); point.back = grid;
    int _d = dir;
    struct { int l, i, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
	POINT_VARIABLES;
	if (condition) {
	  _push (point.level, point.i, 4);
	  if (point.level < point.depth) {
	    _push (point.level + 1, _d == left ? _LEFT : _RIGHT, 0);
	  }
	}
	break;
      }
      case 4: {
        POINT_VARIABLES;
	/* do something */	
@
@def end_foreach_boundary_cell_post()
      }
      }
    }
  }
@

@define foreach_boundary_face(dir) foreach_boundary_cell(dir,true)
@define end_foreach_boundary_face() end_foreach_boundary_cell()

#define update_cache() { if (((Bitree *)grid)->dirty) update_cache_f(); }

@def foreach(clause)     {
  update_cache();
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Bitree point = *((Bitree *)grid); point.back = grid;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < point.leaves.n; _k++) {
    point.i = point.leaves.p[_k].i;
    point.level = point.leaves.p[_k].level;
    POINT_VARIABLES;
@
@define end_foreach() } OMP_END_PARALLEL() }

@def foreach_fine_to_coarse(clause)     {
  update_cache();
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  for (int _l = depth() - 1; _l >= 0; _l--) {
    OMP_PARALLEL()
    Bitree point = *((Bitree *)grid); point.back = grid;
    point.level = _l;
    int _k;
    OMP(omp for schedule(static) clause)
    for (_k = 0; _k < point.active[_l].n; _k++) {
      point.i = point.active[_l].p[_k].i;
      POINT_VARIABLES;
      if (!is_leaf (cell)) {
@
@define end_foreach_fine_to_coarse() } } OMP_END_PARALLEL() } }

@def foreach_level(l) {
  update_cache();
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  int _l = l;
  OMP_PARALLEL()
  Bitree point = *((Bitree *)grid); point.back = grid;
  point.level = _l;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 0; _k < point.active[_l].n; _k++) {
    point.i = point.active[_l].p[_k].i;
    POINT_VARIABLES;
@
@define end_foreach_level() } OMP_END_PARALLEL() }

@def foreach_level_or_leaf(l) {
  for (int _l1 = l; _l1 >= 0; _l1--)
    foreach_level(_l1)
      if (_l1 == l || is_leaf (cell)) {
@
@define end_foreach_level_or_leaf() } end_foreach_level(); }

@define foreach_leaf()            foreach_cell() if (is_leaf (cell)) {
@define end_foreach_leaf()        continue; } end_foreach_cell()

@def foreach_child() {
  int _i = 2*point.i - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2;
  point.level--;
}
@

#if TRASH
# undef trash
# define trash bitree_trash
#endif

void bitree_trash (void * alist)
{
  scalar * list = alist;
  Bitree * q = grid;
#if DYNAMIC
  for (int l = 0; l <= q->depth; l++)
    for (int i = 0; i < _size(l); i++)
      if (q->m[l][i])
	for (scalar s in list)
	  ((double *)(q->m[l][i] + sizeof(Cell)))[s] = undefined;
#else // !DYNAMIC
  int cellsize = sizeof(Cell) + datasize;;
  for (int l = 0; l <= q->depth; l++) {
    char * data = q->m[l] + sizeof(Cell);
    for (int i = 0; i < _size(l); i++, data += cellsize)
      for (scalar s in list)
	((double *)data)[s] = undefined;
  }
#endif // !DYNAMIC
}

#if !DYNAMIC
char * alloc_cells (int l)
{
  int len = _size(l), cellsize = sizeof(Cell) + datasize;
  char * m = calloc (len, cellsize);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  char * data = m + sizeof(Cell);
  int nv = datasize/sizeof(double);
  for (int i = 0; i < len; i++, data += cellsize)
    for (int j = 0; j < nv; j++)
      ((double *)data)[j] = undefined;
  return m;
}
#endif // !DYNAMIC

void alloc_layer (Bitree * p)
{
  Bitree * q = p->back;
  q->depth++; p->depth++;
  q->m = &(q->m[-1]);
#if DYNAMIC
  q->m = realloc(q->m, sizeof (char **)*(q->depth + 2)); 
  q->m = &(q->m[1]); p->m = q->m;
  q->m[q->depth] = calloc (_size(q->depth), sizeof (char *));
#else // !DYNAMIC
  q->m = realloc(q->m, sizeof (char *)*(q->depth + 2)); 
  q->m = &(q->m[1]); p->m = q->m;
  q->m[q->depth] = alloc_cells (q->depth);
#endif
  q->active = realloc (q->active, (q->depth + 1)*sizeof (CacheLevel));
  cache_level_init (&q->active[q->depth]);
  q->halo = realloc (q->halo, (q->depth + 1)*sizeof (CacheLevel));
  cache_level_init (&q->halo[q->depth]);
}

void alloc_children (Bitree * p)
{
  p->back->dirty = true;
  if (p->level == p->depth) alloc_layer(p);
#if DYNAMIC
  char ** m = ((char ***)p->m)[p->level+1];
  Point point = *((Point *)p);
  int l = 0;
  for (int k = - GHOSTS; k < 2 + GHOSTS; k++)
    if (!m[_childindex(k,l)]) {
      m[_childindex(k,l)] = calloc (1, sizeof(Cell) + datasize);
#if TRASH
      char * data = m[_childindex(k,l)] + sizeof(Cell);
      int nv = datasize/sizeof(double);
      for (int j = 0; j < nv; j++)
	((double *)data)[j] = undefined;
#endif
    }
#endif
}

void free_children (Bitree * p)
{
  p->back->dirty = true;
#if DYNAMIC
  char ** m = ((char ***)p->m)[p->level+1];
  Point point = *((Point *)p);
  int l = 0;
  for (int k = - GHOSTS; k < 2 + GHOSTS; k++)
    if (!((Cell *) m[_childindex(k,l)])->neighbors) {
      free (m[_childindex(k,l)]);
      m[_childindex(k,l)] = NULL;
    }
#endif
}

void realloc_scalar (void)
{
  Bitree * q = grid;
#if DYNAMIC
  for (int l = 0; l <= q->depth; l++) {
    size_t len = _size(l);
    for (int i = 0; i < len; i++)
      if (q->m[l][i])
	q->m[l][i] = realloc (q->m[l][i], sizeof(Cell) + datasize);
  }
#else // !DYNAMIC
  size_t oldatasize = sizeof(Cell) + datasize - sizeof(double);
  for (int l = 0; l <= q->depth; l++) {
    size_t len = _size(l);
    q->m[l] = realloc (q->m[l], len*(sizeof(Cell) + datasize));
    char * data = q->m[l] + (len - 1)*oldatasize;
    for (int i = len - 1; i > 0; i--, data -= oldatasize)
      memmove (data + i*sizeof(double), data, oldatasize);
  }
#endif
}

static void update_cache_f (void)
{
  Bitree * q = grid;

  /* empty caches */
  q->leaves.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->halo[l].n = q->active[l].n = 0;

  foreach_cell() {
    if (!is_active (cell)) {
      assert (cell.neighbors > 0);
      /* update halo cache (prolongation) */
      cache_level_append (&q->halo[level], point);
    }
    else {
      if (is_leaf (cell))
	cache_append (&q->leaves, point);
      else if (cell.neighbors > 0)
	/* update halo cache (restriction) */
	cache_level_append (&q->halo[level], point);
      /* update active cache */
      cache_level_append (&q->active[level], point);
    }
  }

  /* update ghost cell flags */
  for (int d = 0; d < nboundary; d++)
    foreach_boundary_cell (d, true) {
      neighbor(ghost).flags = fghost;
      if (!is_active (cell))
	continue;
    }

  q->dirty = false;
}

@def foreach_halo_level(_l) {
  update_cache();
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Bitree point = *((Bitree *)grid); point.back = grid;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 0; _k < point.halo[_l].n; _k++) {
    point.i = point.halo[_l].p[_k].i;
    point.level = _l;
    POINT_VARIABLES;
@
@define end_foreach_halo_level() } OMP_END_PARALLEL() }

@def foreach_halo_levels(start,cond,inc) {
  for (int _l = start; _l cond; _l inc)
    foreach_halo_level(_l)
@
@define end_foreach_halo_levels() end_foreach_halo_level() }

/* breadth-first traversal of halos from coarse to fine */
@def foreach_halo_coarse_to_fine(depth1) {
  foreach_halo_levels (0, <= depth1, ++)
    if (!is_active(cell)) {
@
@define end_foreach_halo_coarse_to_fine() } end_foreach_halo_levels() }

/* breadth-first traversal of halos from fine to coarse */
@def foreach_halo_fine_to_coarse()
  foreach_halo_levels (depth() - 1, >= 0, --)
    if (is_active(cell)) {
@
@define end_foreach_halo_fine_to_coarse() } end_foreach_halo_levels()

@def foreach_halo_vertex()
  foreach_halo_levels(0, <= depth(), ++)
  if (allocated(-1,0) && is_leaf(neighbor(-1,0))) {
@
@define end_foreach_halo_vertex() } end_foreach_halo_levels()

Point refine_cell (Point point, scalar * list);

void free_grid (void)
{
  if (!grid)
    return;
  Bitree * q = grid;
  free (q->leaves.p);
  for (int l = 0; l <= q->depth; l++) {
#if DYNAMIC
    for (int i = 0; i < _size(l); i++)
      free (q->m[l][i]);
#endif
    free (q->m[l]);
    free (q->halo[l].p);
    free (q->active[l].p);
  }
  q->m = &(q->m[-1]);
  free (q->m);
  free (q->halo);
  free (q->active);
  free (q);
  grid = NULL;
}

void init_grid (int n)
{
  Bitree * q = grid;
  if (q && n == 1 << q->depth)
    return;
  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "bitree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  q = malloc (sizeof (Bitree));
  q->depth = 0; q->i = q->j = GHOSTS; q->level = 0.;
  q->m = malloc(sizeof (char *)*2);
  /* make sure we don't try to access level -1 */
  q->m[0] = NULL; q->m = &(q->m[1]);
  /* initialise the root cell */
#if DYNAMIC
  int len = _size(0);
  q->m[0] = calloc (len, sizeof (char *));
  for (int i = 0; i < len; i++)
    q->m[0][i] = calloc (1, sizeof(Cell) + datasize);
#else
  q->m[0] = alloc_cells (0);
#endif
  CELL(q->m, 0, GHOSTS).flags |= (leaf | active);
  CELL(q->m, 0, GHOSTS).neighbors = 1; // only itself as neighbor
  cache_init (&q->leaves);
  q->active = calloc (1, sizeof (CacheLevel));
  q->halo = calloc (1, sizeof (CacheLevel));
  q->dirty = true;
  grid = q;
  N = 1 << depth;
  while (depth--)
    foreach_leaf()
      point = refine_cell (point, NULL);
  update_cache();
  trash (all);
  init_events();
}

void output_cells (FILE * fp);

#include "quadtree-common.h"

void bitree_methods()
{
  quadtree_methods();
}
