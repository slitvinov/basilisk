#include "mempool.h"

#define GRIDNAME "Quadtree"
#define dimension 2

#define TWO_ONE 1 // enforce 2:1 refinement ratio
#define GHOSTS  2

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

typedef struct {
  short flags;
  short neighbors; // number of refined neighbors in a 3x3 neighborhood
  int pid;         // process id
} Cell;

enum {
  active  = 1 << 0,
  leaf    = 1 << 1,
  refined = 1 << 2,
  halo    = 1 << 3,
  border  = 1 << 4,
  remote_leaf = 1 << 5,
  fboundary = 1 << 6,
  user    = 7,

  face_x = 1 << 0,
  face_y = 1 << 1
};

@define is_active(cell)  ((cell).flags & active)
@define is_leaf(cell)    ((cell).flags & leaf)
@define is_coarse()      ((cell).neighbors > 0)
@define is_border(cell)  ((cell).flags & border)
@define is_local(cell)   ((cell).pid == pid())

@if _MPI
@ define is_remote_leaf(cell) ((cell).flags & remote_leaf)
@else
@ define is_remote_leaf(cell) false
@endif

// Caches

typedef struct {
  int i, j;
} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i, j;
  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;

// Layer

typedef struct {
  char *** m;     // the 2D array of data
  Mempool * pool; // the memory pool actually holding the data
  int * nr;   // the number of allocated rows for each column
  int nc;     // the number of allocated columns
  int len;    // the (1D) size of the array
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*GHOSTS;
}

static size_t poolsize (size_t depth, size_t size)
{
  // the maximum amount of data at a given level
  return sq(_size(depth))*size;
}

static Layer * new_layer (int depth)
{
  Layer * l = malloc (sizeof (Layer));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL; // the root layer does not use a pool
  else {
    size_t size = sizeof(Cell) + datasize;
    // the block size is 4*size because we allocate 4 children at a time
    l->pool = mempool_new (poolsize (depth, size), 4*size);
  }
  l->m = calloc (l->len, sizeof (char *));
  l->nr = calloc (l->len, sizeof (int));
  l->nc = 0;
  return l;
}

static void layer_add_row (Layer * l, int i)
{
  if (!l->m[i]) {
    l->m[i] = calloc (l->len, sizeof (char *));
    l->nc++;
  }
  l->nr[i]++;
}

static bool layer_remove_row (Layer * l, int i)
{
  if (--l->nr[i] == 0) {
    free (l->m[i]);
    l->m[i] = NULL;
    if (--l->nc == 0) {
      // fixme: need global depth in parallel
      return false;

      free (l->m);
      free (l->nr);
      free (l);
      return true; // layer has been destroyed
    }
  }
  return false;
}

// Quadtree

typedef struct {
  Layer ** L; /* the grids at each level */
  int depth;  /* the maximum depth of the tree */

  Cache        leaves;   /* leaf indices */
  Cache        faces;    /* face indices */
  Cache        vertices; /* vertex indices */
  Cache        refined;  /* refined cells */
  CacheLevel * active;   /* active cells indices for each level */
  CacheLevel * prolongation; /* halo prolongation indices for each level */
  CacheLevel * restriction;  /* halo restriction indices for each level */
  CacheLevel * boundary;  /* boundary indices for each level */
  CacheLevel   coarsened; /* coarsened cells */
  
  bool dirty;       /* whether caches should be updated */
} Quadtree;

#define quadtree ((Quadtree *)grid)

struct _Point {
  int i, j, level;  /* the current cell index and level */
};
static Point last_point;

static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (IndexLevel)*c->nm);
  }
  c->p[c->n].i = p.i;
  c->p[c->n].j = p.j;
  c->n++;
}

static void cache_append (Cache * c, Point p, int k, int l, short flags)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (Index)*c->nm);
  }
  c->p[c->n].i = p.i + k;
  c->p[c->n].j = p.j + l;
  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

@define cache_append_x(c,p,k,l,flags) cache_append(c,p,k,l,flags)
@define cache_append_y(c,p,k,l,flags) cache_append(c,p,l,k,flags)

/* low-level memory management */
@def allocated(k,l,n) (quadtree->L[point.level]->m[point.i+k] &&
		       quadtree->L[point.level]->m[point.i+k][point.j+l])
@
@define NEIGHBOR(k,l)	(quadtree->L[point.level]->m[point.i+k][point.j+l])
@def PARENT(k,l) (quadtree->L[point.level-1]->m[(point.i+GHOSTS)/2+k]
		  [(point.j+GHOSTS)/2+l])
@
@def allocated_child(k,l,n)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k]
   && quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k][2*point.j-GHOSTS+l])
@			   
@def CHILD(k,l)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k]
		 [2*point.j-GHOSTS+l])
@
@define CELL(m) (*((Cell *)(m)))

/***** Multigrid macros *****/
@define depth()        (quadtree->depth)
@define aparent(k,l,n) CELL(PARENT(k,l))
@define child(k,l,n)   CELL(CHILD(k,l))

/***** Quadtree macros ****/
@define NN              (1 << point.level)
@define cell		CELL(NEIGHBOR(0,0))
@define neighbor(k,l,n)	CELL(NEIGHBOR(k,l))

/***** Data macros *****/
@define data(k,l,n)     ((double *) (NEIGHBOR(k,l) + sizeof(Cell)))
@define fine(a,k,l,n)   ((double *) (CHILD(k,l) + sizeof(Cell)))[a]
@define coarse(a,k,l,n) ((double *) (PARENT(k,l) + sizeof(Cell)))[a]

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
  if (point.level == quadtree->depth) {
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

@def foreach_cell()
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */
    if (is_active(cell))
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, point.j, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
	POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell()
        if (point.level < quadtree->depth) {
	  _push (point.level, point.i, point.j, 1);
          _push (point.level + 1, _LEFT, _TOP, 0);
        }
        break;
      }
      case 1: _push (point.level, point.i, point.j, 2);
              _push (point.level + 1, _RIGHT, _TOP,    0); break;
      case 2: _push (point.level, point.i, point.j, 3);
              _push (point.level + 1, _LEFT,  _BOTTOM, 0); break;
      case 3: _push (point.level + 1, _RIGHT, _BOTTOM, 0); break;
      }
    }
  }
@

@def foreach_cell_post(condition)
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, point.j, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
	if (point.level == quadtree->depth) {
	  _push (point.level, point.i, point.j, 4);
	}
	else {
	  _push (point.level, point.i, point.j, 1);
	  if (condition)
	    _push (point.level + 1, _LEFT, _TOP, 0);
	}
	break;
      }
      case 1:
	_push (point.level, point.i, point.j, 2);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP,    0);
	break;
      case 2:
	_push (point.level, point.i, point.j, 3);
	if (condition)
	  _push (point.level + 1, _LEFT,  _BOTTOM, 0);
	break;
      case 3:
	_push (point.level, point.i, point.j, 4);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, 0);
	break;
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

@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 2; _k++)
    for (int _l = 0; _l < 2; _l++) {
      point.i = _i + _k; point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
}
@

@def foreach_child_direction(d) {
  static int _ii[4][2] = {{1,1}, {0,0}, {0,1}, {0,1}};
  static int _jj[4][2] = {{0,1}, {0,1}, {1,1}, {0,0}};
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _ii[d][_k]; point.j = _j + _jj[d][_k];
    POINT_VARIABLES;
@
@def end_foreach_child_direction()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
}
@

#define update_cache() { if (quadtree->dirty) update_cache_f(); }

#define is_prolongation(cell) (!is_leaf(cell) && !cell.neighbors && \
			       cell.pid >= 0)
#define is_boundary(cell) (cell.pid < 0)

@def foreach_cache(_cache,clause) {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Point point = {GHOSTS,GHOSTS,0};
  int _k; short _flags; NOT_UNUSED(_flags);
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
    point.level = _cache.p[_k].level;
    _flags = _cache.p[_k].flags;
    POINT_VARIABLES;
@
@define end_foreach_cache() } OMP_END_PARALLEL() }

@def foreach_cache_level(_cache,_l,clause) {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Point point = {GHOSTS,GHOSTS,0};
  point.level = _l;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
    POINT_VARIABLES;
@
@define end_foreach_cache_level() } OMP_END_PARALLEL() }

@def foreach_boundary(l) {
  update_cache();
  CacheLevel _boundary = quadtree->boundary[l];
  foreach_cache_level (_boundary,l,)
@
@define end_foreach_boundary() end_foreach_cache_level(); }

@def foreach_halo(name,_l) {
  update_cache();
  CacheLevel _cache = quadtree->name[_l];
  foreach_cache_level (_cache, _l,)
@
@define end_foreach_halo() end_foreach_cache_level(); }

static void update_cache_f (void)
{
  Quadtree * q = grid;

  /* empty caches */
  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->restriction[l].n = q->boundary[l].n = 0;

  foreach_cell() {
    if (is_local(cell)) {
      // active cells
      assert (is_active(cell));
      cache_level_append (&q->active[level], point);
    }
    // boundaries
    if (!is_boundary(cell))
      for (int k = -1; k <= 1; k++)
	for (int l = -1; l <= 1; l++)
	  if (is_boundary(neighbor(k,l)) && !(neighbor(k,l).flags & fboundary)) {
	    point.i += k; point.j += l;
	    cache_level_append (&q->boundary[level], point);
	    cell.flags |= fboundary;
	    point.i -= k; point.j -= l;
	  }
    if (is_leaf (cell)) {
      if (is_local(cell)) {
	cache_append (&q->leaves, point, 0, 0, 0);
	// faces
	unsigned flags = 0;
	foreach_dimension()
	  if (is_boundary(neighbor(-1,0)) ||
	      is_prolongation(neighbor(-1,0)) || is_leaf(neighbor(-1,0)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, 0, 0, flags);
	foreach_dimension()
	  if (is_boundary(neighbor(1,0)) || is_prolongation(neighbor(1,0)) ||
	      (!is_local(neighbor(1,0)) && is_leaf(neighbor(1,0))))
	    cache_append_x (&q->faces, point, 1, 0, face_x);
	// vertices
	cache_append (&q->vertices, point, 0, 0, 0);
	if (!is_leaf(neighbor(1,1)) || !is_local(neighbor(1,1)))
	  cache_append (&q->vertices, point, 1, 1, 0);
	foreach_dimension()
	  if ((!is_leaf(neighbor(0,-1)) || !is_local(neighbor(0,-1))) &&
	      (!is_leaf(neighbor(1,0)) || !is_local(neighbor(1,0))))
	    cache_append_x (&q->vertices, point, 1, 0, 0);
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
	cell.flags &= ~halo;
      }
      else if (!is_boundary(cell) || is_local(aparent(0,0))) { // non-local
	// faces
	unsigned flags = 0;
	foreach_dimension()
	  if (allocated(-1,0) &&
	      is_local(neighbor(-1,0)) && is_prolongation(neighbor(-1,0)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, 0, 0, flags);
	foreach_dimension()
	  if (allocated(1,0) &&
	      is_local(neighbor(1,0)) && is_prolongation(neighbor(1,0)))
	    cache_append_x (&q->faces, point, 1, 0, face_x);
	// vertices: fixme: this is too complicated
	foreach_dimension()
	  for (int i = -1; i <= 1; i += 2)
	    if (allocated(i,0) &&
		is_local(neighbor(i,0)) && is_prolongation(neighbor(i,0))) {
	      cache_append_x (&q->vertices, point, (i + 1)/2, 0, 0);
	      cache_append_x (&q->vertices, point, (i + 1)/2, 1, 0);
	    }
      }
      continue;
    }
    else { // not a leaf
      bool restriction =
	level > 0 &&
	(aparent(0,0).flags & halo);
      for (int k = -GHOSTS; k <= GHOSTS && !restriction; k++)
	for (int l = -GHOSTS; l <= GHOSTS && !restriction; l++)
	  if (allocated(k,l) &&
	      ((is_local(neighbor(k,l)) && is_leaf(neighbor(k,l))) ||
	       is_remote_leaf(neighbor(k,l))))
	    restriction = true;
      if (restriction) {
	// halo restriction
	cell.flags |= halo;
	if (is_local(cell))
	  cache_level_append (&q->restriction[level], point);
      }
      else
	cell.flags &= ~halo;
    }
  }
  
  q->dirty = false;

  for (int l = depth(); l >= 0; l--) {
    foreach_boundary (l)
      cell.flags &= ~(fboundary|halo);
    // we mark boundary cells which are required for halo restriction
    foreach_halo (restriction, l)
      for (int k = 0; k <= 2; k++)
	for (int l = 0; l <= 2; l++)
	  if (is_boundary(child(k,l)))
	    child(k,l).flags |= halo;
  }
}

@define foreach(clause) update_cache(); foreach_cache(quadtree->leaves, clause)
@define end_foreach()   end_foreach_cache()

@def foreach_face_generic(clause)
  update_cache();
  foreach_cache(quadtree->faces, clause) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() (_flags & face_x)
@define is_face_y() (_flags & face_y)

@def foreach_vertex(clause)
  update_cache();
  foreach_cache(quadtree->vertices, clause) {
    x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_cache()

@def foreach_fine_to_coarse(clause)     {
  update_cache();
  for (int _l = depth() - 1; _l >= 0; _l--) {
    CacheLevel _active = quadtree->active[_l];
    foreach_cache_level (_active,_l,clause)
      if (!is_leaf (cell)) {
@
@define end_foreach_fine_to_coarse() } end_foreach_cache_level(); } }

@def foreach_level(l) {
  update_cache();
  CacheLevel _active = quadtree->active[l];
  foreach_cache_level (_active,l,)
@
@define end_foreach_level() end_foreach_cache_level(); }

@define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
@define end_foreach_coarse_level() } end_foreach_level()

@def foreach_level_or_leaf(l) {
  for (int _l1 = l; _l1 >= 0; _l1--)
    foreach_level(_l1)
      if (_l1 == l || is_leaf (cell)) {
@
@define end_foreach_level_or_leaf() } end_foreach_level(); }

@def foreach_leaf() foreach_cell()
  if (is_leaf (cell)) {
    if (is_active(cell)) {
@
@define end_foreach_leaf() } continue; } end_foreach_cell()

@if TRASH
@ undef trash
@ define trash quadtree_trash
@endif

void quadtree_trash (void * alist)
{
  scalar * list = alist;
  Quadtree * q = grid;
  /* low-level memory management */
  for (int l = 0; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])
	for (int j = 0; j < L->len; j++)
	  if (L->m[i][j])
	    for (scalar s in list)
	      ((double *)(L->m[i][j] + sizeof(Cell)))[s] = undefined;
  }
}

@def cache_level_resize(name, a)
{
  for (int i = 0; i <= q->depth - a; i++)
    free (q->name[i].p);
  free (q->name);
  q->name = calloc (q->depth + 1, sizeof (CacheLevel));
}
@

static void alloc_layer (void)
{
  // fixme: need global depth in parallel

  Quadtree * q = quadtree;
  q->depth++;
  /* low-level memory management */
  q->L = &(q->L[-1]);
  q->L = realloc(q->L, sizeof (Layer *)*(q->depth + 2));
  q->L = &(q->L[1]);
  q->L[q->depth] = new_layer (q->depth);
  cache_level_resize (active, +1);
  cache_level_resize (prolongation, +1);
  cache_level_resize (restriction, +1);
  cache_level_resize (boundary, +1);
}

static void alloc_children (Point point, int i, int j)
{
  if (point.level == quadtree->depth)
    alloc_layer();
  point.i += i; point.j += j;

  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = mempool_alloc0 (L->pool);
  for (int k = 0; k < 2; k++) {
    layer_add_row (L, 2*point.i - GHOSTS + k);
    for (int l = 0; l < 2; l++) {
      assert (!CHILD(k,l));
      CHILD(k,l) = b;
      b += len;
    }
  }

  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      child(k,l).pid = cell.pid;

@if TRASH
  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      for (scalar s in all)
	fine(s,k,l) = undefined;
@endif
}

void increment_neighbors (Point point)
{
  quadtree->dirty = true;
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++) {
      if (neighbor(k,l).neighbors == 0)
	alloc_children (point, k, l);
      neighbor(k,l).neighbors++;
    }
}

static void free_layer (void)
{
  Quadtree * q = quadtree;
  q->depth--;
  /* low-level memory management */
  q->L = &(q->L[-1]);
  q->L = realloc(q->L, sizeof (Layer *)*(q->depth + 2));
  q->L = &(q->L[1]);
  cache_level_resize (active, -1);
  cache_level_resize (prolongation, -1);
  cache_level_resize (restriction, -1);
  cache_level_resize (boundary, -1);
}

static void free_children (Point point, int i, int j)
{
  point.i += i; point.j += j;
  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  assert (CHILD(0,0));
  mempool_free (L->pool, CHILD(0,0));
  for (int k = 0; k < 2; k++) {
    for (int l = 0; l < 2; l++)
      CHILD(k,l) = NULL;
    if (layer_remove_row (L, 2*point.i - GHOSTS + k)) {
      assert (point.level + 1 == quadtree->depth);
      free_layer();
    }
  }
}

void decrement_neighbors (Point point)
{
  quadtree->dirty = true;
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
      if (allocated(k,l)) {
	neighbor(k,l).neighbors--;
	if (neighbor(k,l).neighbors == 0)
	  free_children (point,k,l);
      }
}

void realloc_scalar (void)
{
  /* low-level memory management */
  Quadtree * q = grid;
  size_t newlen = sizeof(Cell) + datasize;
  size_t oldlen = newlen - sizeof(double);
  /* the root level is allocated differently */
  size_t len = _size(0);
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
      q->L[0]->m[i][j] = realloc (q->L[0]->m[i][j], newlen);
  /* all other levels */
  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    char *** m = L->m;
    size_t len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), 4*newlen);
    for (int i = 0; i < len; i += 2)
      if (m[i])
	for (int j = 0; j < len; j += 2)
	  if (m[i][j]) {
	    char * new = mempool_alloc (L->pool);
	    for (int k = 0; k < 2; k++)
	      for (int o = 0; o < 2; o++) {
		memcpy (new, m[i+k][j+o], oldlen);
		m[i+k][j+o] = new;
		new += newlen;
	      }
	  }
    mempool_destroy (oldpool);
  }
}

/* Boundaries */

@def foreach_normal_neighbor(cond)
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if ((i == 0 || j == 0) && allocated(i,j) && !is_boundary(neighbor(i,j))) {
	Point neighbor = {point.i + i, point.j + j, point.level};
	if (cond (neighbor)) {
@
@def end_foreach_normal_neighbor()
        }
      }
@

@def foreach_normal_neighbor_x(cond)
  for (int i = -1; i <= 1; i += 2)
    if (allocated(i,0) && !is_boundary(neighbor(i,0))) {
      Point neighbor = {point.i + i, point.j, point.level};
      if (cond (neighbor)) {
@
@def end_foreach_normal_neighbor_x()
      }
    }
@

@def foreach_normal_neighbor_y(cond)
  for (int i = -1; i <= 1; i += 2)
    if (allocated(0,i) && !is_boundary(neighbor(0,i))) {
      Point neighbor = {point.i, point.j + i, point.level};
      if (cond (neighbor)) {
@
@def end_foreach_normal_neighbor_y()
      }
    }
@

@def foreach_tangential_neighbor_x(cond)
  for (int j = -1; j <= 1; j += 2) {
    Point neighbor = {point.i, point.j + j, point.level};
    Point n2 = {point.i - 1, point.j + j, point.level};
    if ((allocated(0,j) && !is_boundary(neighbor(0,j)) && cond(neighbor)) ||
	(allocated(-1,j) && !is_boundary(neighbor(-1,j)) && cond(n2))) {
@
@def end_foreach_tangential_neighbor_x()
    }
  }
@

@def foreach_tangential_neighbor_y(cond)
  for (int j = -1; j <= 1; j += 2) {
    Point neighbor = {point.i + j, point.j, point.level};
    Point n2 = {point.i + j, point.j - 1, point.level};
    if ((allocated(j,0) && !is_boundary(neighbor(j,0)) && cond(neighbor)) ||
	(allocated(j,-1) && !is_boundary(neighbor(j,-1)) && cond(n2))) {
@
@def end_foreach_tangential_neighbor_y()
    }
  }
@

#define bid(cell) (- cell.pid - 1)

static int boundary_scalar (Point point, bool (*cond)(Point), scalar * list)
{
  if (!list)
    return 1;
  int nc = 0, id = bid(cell);
  for (scalar s in list)
    s[] = 0.;
  /* look first in normal directions */
  foreach_normal_neighbor (cond) {
    for (scalar s in list)
      s[] += s.boundary[id] (neighbor, point, s);
    nc++;
  }
  if (nc > 0)
    for (scalar s in list)
      s[] /= nc;
  return nc;
}

@define VN _attribute[s].v.x
@define VT _attribute[s].v.y
  
foreach_dimension()
static int boundary_vector_x (Point point, bool (*cond)(Point), scalar * list)
{
  if (!list)
    return 1;
  int nc = 0, id = bid(cell);
  for (scalar s in list)
    s[] = 0.;
  /* look first in normal directions */
  foreach_normal_neighbor_x (cond) {
    for (scalar s in list) {
      scalar n = VN;
      s[] += n.boundary[id] (neighbor, point, s);
    }
    nc++;
  }
  if (!nc)
    /* look in tangential directions */
    foreach_normal_neighbor_y (cond) {
      for (scalar s in list) {
	scalar t = VT;
	s[] += t.boundary[id] (neighbor, point, s);
      }
      nc++;
    }
  if (nc > 0)
    for (scalar s in list)
      s[] /= nc;
  return nc;
}

static void boundary_normal_x (Point point, bool (*cond) (Point), scalar * list)
{
  if (!list)
    return;
  int id = bid(cell);
  foreach_normal_neighbor_x (cond)
    for (scalar s in list) {
      scalar n = VN;
      if (s.normal && n.boundary[id])
	s[(i + 1)/2,0] = n.boundary[id] (neighbor, point, s);
    }
}

static void boundary_normal_y (Point point, bool (*cond) (Point), scalar * list)
{
  if (!list)
    return;
  int id = bid(cell);
  foreach_normal_neighbor_y (cond)
    for (scalar s in list) {
      scalar n = VN;
      if (s.normal && n.boundary[id])
	s[0,(i + 1)/2] = n.boundary[id] (neighbor, point, s);
    }
}

foreach_dimension()
static void boundary_tangential_x (Point point, bool (*cond)(Point),
				   scalar * list)
{
  if (!list)
    return;
  if (allocated(-1,0) && neighbor(-1,0).pid < 0) {
    int nc = 0, id = bid(cell);
    for (scalar s in list)
      s[] = 0.;
    foreach_tangential_neighbor_x (cond) {
      for (scalar s in list) {
	scalar t = VT;
	s[] += t.boundary[id] (neighbor, point, s);
      }
      nc++;
    }
    if (nc > 0) {
      for (scalar s in list)
	s[] /= nc;
    }
    else
      for (scalar s in list)
	s[] = 0.; // fixme: this should be undefined
  }
}  

static void boundary_diagonal (Point point, bool (*cond)(Point), scalar * list)
{
  if (!list)
    return;
  int nc = 0;
  for (int k = -1; k <= 1; k += 2)
    for (int l = -1; l <= 1; l += 2)
      if (allocated(k,l) && neighbor(k,l).pid >= 0) {
	Point neighbor = {point.i + k, point.j + l, point.level};
	if (cond (neighbor)) {
	  for (scalar s in list)
	    s[] += s[k,0] + s[0,l] - s[k,l];
	  nc++;
	}
      }
  if (nc > 0) {
    for (scalar s in list)
      s[] /= nc;
  }
  else
    for (scalar s in list)
      s[] = undefined;
}
 
#define box_boundaries(l, cond1, cond, list)				\
do {									\
  scalar * lists = NULL, * list_x = NULL, * list_y = NULL;		\
  scalar * listf_x = NULL, * listf_y = NULL;				\
  for (scalar s in list)						\
    if (!is_constant(s)) {						\
      if (s.v.x < 0)							\
	lists = list_add (lists, s);					\
      else {								\
	if (s.v.x == s) {						\
	  if (s.face)							\
	    listf_x = list_add (listf_x, s);				\
	  else								\
	    list_x = list_add (list_x, s);				\
	}								\
	else if (s.face)						\
	  listf_y = list_add (listf_y, s);				\
	else								\
	  list_y = list_add (list_y, s);				\
      }									\
    }									\
  /* normal/tangential directions */					\
  int diagonal = 1 << user, diagonal_x = 2*diagonal, diagonal_y = 4*diagonal; \
  foreach_boundary(l)							\
    if (cond1) {							\
      if (!boundary_scalar (point, cond, lists))			\
	cell.flags |= diagonal;						\
      foreach_dimension() {						\
	if (!boundary_vector_x (point, cond, list_x))			\
	  cell.flags |= diagonal_x;					\
	boundary_normal_x (point, cond, listf_x);			\
      }									\
    }									\
  /* diagonal and tangential (face) directions */			\
  foreach_boundary(l)							\
    if (cond1) {							\
      if (cell.flags & diagonal) {					\
	boundary_diagonal (point, cond, lists);				\
	cell.flags &= ~diagonal;					\
      }									\
      foreach_dimension() {						\
	if (cell.flags & diagonal_x) {					\
	  boundary_diagonal (point, cond, list_x);			\
	  cell.flags &= ~diagonal_x;					\
	}								\
	boundary_tangential_x (point, cond, listf_x);			\
      }									\
    }									\
  free (lists);								\
  foreach_dimension() {							\
    free (list_x);							\
    free (listf_x);							\
  }									\
 } while(0)

static bool retrue (Point point) { return true; }

static bool retleaf (Point point) { return is_leaf(cell); }

static bool retleafhalo (Point point) {
  // leaf or prolongation or restriction halo
  return is_leaf(cell) || !cell.neighbors || (cell.flags & halo);
}
 
static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  /* we disable floating-point-exceptions to avoid having to deal with
     undefined operations in non-trivial boundary conditions. */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  if (l < 0)
    for (l = 0; l <= depth(); l++)
      box_boundaries (l, true, retleaf, list);
  else
    box_boundaries (l, true, retrue, list);
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

static void box_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  if (l == 0)
    return;
  /* we disable floating-point-exceptions to avoid having to deal with
     undefined operations in non-trivial boundary conditions. */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  box_boundaries (l, cell.flags & halo, retrue, list);
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

static void box_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, 
					    int l, int depth)
{
  /* we disable floating-point-exceptions to avoid having to deal with
     undefined operations in non-trivial boundary conditions. */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  if (l == depth)
    box_boundaries (l, true, retrue, list);
  else
    box_boundaries (l, true, retleafhalo, list);
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
 
#define mask(func) {					\
  foreach_cell_post(!is_leaf(cell)) {			\
    if (is_leaf(cell)) {				\
      int bid = (func);					\
      if (bid >= 0)					\
	cell.pid = - bid - 1;				\
    }							\
    else { /* not a leaf */				\
      int pid = -1;					\
      foreach_child() {					\
	if (cell.pid >= 0 || pid < 0)			\
	  pid = cell.pid;				\
      }							\
      cell.pid = pid;					\
      if (pid < 0) {					\
	/* fixme: call coarsen_cell()? */		\
	cell.flags |= leaf;				\
	decrement_neighbors (point);			\
	if (cell.neighbors)				\
	  for (int k = 0; k < 2; k++)			\
	    for (int l = 0; l < 2; l++) {		\
	      child(k,l).flags &= ~(leaf|active);	\
	      child(k,l).pid = cell.pid;		\
	    }						\
      }							\
    }							\
  }							\
  quadtree->dirty = true;				\
}

/* Periodic boundaries */

static double periodic_bc (Point, Point, scalar);
  
foreach_dimension()
static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face) {
	scalar vt = VT;
	if (vt.boundary[right] == periodic_bc)
	  list1 = list_add (list1, s);
      }
      else if (s.boundary[right] == periodic_bc)
	list1 = list_add (list1, s);
    }
  if (!list1)
    return;

  OMP_PARALLEL();
  Point point = {0,0};
  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;
  int j;
  OMP(omp for schedule(static))
  for (j = 0; j < n + 2*GHOSTS; j++) {
    for (int i = 0; i < GHOSTS; i++)
      if (allocated(i + n,j))
	for (scalar s in list1)
	  s[i,j] = s[i + n,j];
    for (int i = n + GHOSTS; i < n + 2*GHOSTS; i++)
      if (allocated(i - n,j))
	for (scalar s in list1)
	  s[i,j] = s[i - n,j];
  }
  OMP_END_PARALLEL();
  
  free (list1);
}

@undef VN
@undef VT

foreach_dimension()
static void periodic_boundary_halo_prolongation_x (const Boundary * b,
						   scalar * list, 
						   int l, int depth)
{
  periodic_boundary_level_x (b, list, l);
}

void refine_cell (Point point, scalar * list, int flag, Cache * refined);

static void free_cache (CacheLevel * c)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++)
    free (c[l].p);
  free (c);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Quadtree * q = grid;
  free (q->leaves.p);
  free (q->faces.p);
  free (q->vertices.p);
  free (q->refined.p);
  free (q->coarsened.p);
  /* low-level memory management */
  /* the root level is allocated differently */
  Layer * L = q->L[0];
  for (int i = 0; i < L->len; i++) {
    for (int j = 0; j < L->len; j++)
      free (L->m[i][j]);
    free (L->m[i]);
  }
  free (L->m);
  free (L->nr);
  free (L);
  /* all other levels */
  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])
	free (L->m[i]);
    mempool_destroy (L->pool);
    free (L->m);
    free (L->nr);
    free (L);
  }
  q->L = &(q->L[-1]);
  free (q->L);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->restriction);
  free_cache (q->boundary);
  free (q);
  grid = NULL;
}

trace
void init_grid (int n)
{
  // check 64 bits structure alignment
  assert (sizeof(Cell) % 8 == 0);
  
  Quadtree * q = grid;
  if (q && n == 1 << q->depth)
    return;
  free_grid();

  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  q = calloc (1, sizeof (Quadtree));
  q->depth = 0;

  /* low-level memory management */
  q->L = malloc(sizeof (Layer *)*2);
  /* make sure we don't try to access level -1 */
  q->L[0] = NULL; q->L = &(q->L[1]);
  /* initialise the root cell */
  Layer * L = new_layer (0);
  q->L[0] = L;
  for (int i = 0; i < L->len; i++) {
    layer_add_row (L, i);
    for (int j = 0; j < L->len; j++)
      L->m[i][j] = calloc (1, sizeof(Cell) + datasize);
  }
  CELL(L->m[GHOSTS][GHOSTS]).flags |= (leaf|active);
  for (int k = -GHOSTS; k <= GHOSTS; k++)
    for (int l = -GHOSTS; l <= GHOSTS; l++)
      CELL(L->m[GHOSTS+k][GHOSTS+l]).pid =
	(k < 0 ? -1 - left :
	 k > 0 ? -1 - right :
	 l > 0 ? -1 - top :
	 l < 0 ? -1 - bottom :
	 pid());
  q->active = calloc (1, sizeof (CacheLevel));
  q->prolongation = calloc (1, sizeof (CacheLevel));
  q->restriction = calloc (1, sizeof (CacheLevel));
  q->boundary = calloc (1, sizeof (CacheLevel));
  q->dirty = true;
  grid = q;
  N = 1 << depth;
  while (depth--)
    foreach_leaf()
      refine_cell (point, NULL, 0, NULL);
@if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
@endif
  update_cache();
  trash (all);
  // boundaries
  Boundary * b = calloc (1, sizeof (Boundary));
  b->level = b->restriction = box_boundary_level;
  b->halo_restriction  = box_boundary_halo_restriction;
  b->halo_prolongation = box_boundary_halo_prolongation;
  add_boundary (b);
  // periodic boundaries
  foreach_dimension() {
    Boundary * b = calloc (1, sizeof (Boundary));
    b->level = b->restriction = periodic_boundary_level_x;
    b->halo_prolongation = periodic_boundary_halo_prolongation_x;
    add_boundary (b);
  }
@if _MPI
  void mpi_partitioning();
  if (N > 1)
    mpi_partitioning();
@endif
  init_events();
}

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

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + GHOSTS;
    point.j = (p.y - Y0)/L0*n + GHOSTS;
    if (point.i >= 0 && point.i < n + 2*GHOSTS) {
      if (allocated(0,0) && is_local(cell) && is_leaf(cell))
	return point;
    }
    else
      break;
  }
  Point point = { .level = -1 };
  return point;
}

#include "quadtree-common.h"

@if _MPI
#include "quadtree-mpi.h"
@endif
