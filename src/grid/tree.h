#include "mempool.h"

#define TWO_ONE 1 // enforce 2:1 refinement ratio
#define GHOSTS  2

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

#if dimension == 3
# define K (point.k - GHOSTS)
#endif

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
#if dimension == 3
  , face_z = 1 << 2
#endif
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
#if dimension == 3
  int k;
#endif
} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i, j;
#if dimension == 3
  int k;
#endif  
  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;

// Refcounted array

static void * new_refarray (size_t len, size_t size) {
  return calloc (len + 1, size);
}

static void refarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)++;
}

static bool unrefarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)--;
  if (*refcount == 0) {
    free (p);
    return true;
  }
  return false;
}

// Layer

typedef struct {
#if dimension == 2
  char *** m;     // the 2D array of data
#else // dimension == 3
  char **** m;    // the 3D array of data
#endif
  Mempool * pool; // the memory pool actually holding the data
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
#if dimension == 2
  return sq(_size(depth))*size;
#else
  return cube(_size(depth))*size;
#endif
}

static Layer * new_layer (int depth)
{
  Layer * l = malloc (sizeof (Layer));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL; // the root layer does not use a pool
  else {
    size_t size = sizeof(Cell) + datasize;
    // the block size is 2^dimension*size because we allocate
    // 2^dimension children at a time
    l->pool = mempool_new (poolsize (depth, size), (1 << dimension)*size);
  }
  l->m = calloc (l->len, sizeof (char *));
  l->nc = 0;
  return l;
}

#if dimension == 2
static void layer_add_row (Layer * l, int i)
{
  if (!l->m[i]) {
    l->m[i] = new_refarray (l->len, sizeof (char *));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));
}

static bool layer_remove_row (Layer * l, int i)
{
  if (unrefarray (l->m[i], l->len, sizeof (char *))) {
    l->m[i] = NULL;
    if (--l->nc == 0) {
      // fixme: need global depth in parallel
      return false;

      free (l->m);
      free (l);
      return true; // layer has been destroyed
    }
  }
  return false;
}
#else // dimension == 3
static void layer_add_row (Layer * l, int i, int j)
{
  if (!l->m[i]) {
    l->m[i] = new_refarray (l->len, sizeof (char *));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));
  if (!l->m[i][j])
    l->m[i][j] = new_refarray (l->len, sizeof (char *));
  refarray (l->m[i][j], l->len, sizeof(char *));  
}

static bool layer_remove_row (Layer * l, int i, int j)
{
  if (unrefarray (l->m[i][j], l->len, sizeof (char *))) {
    l->m[i][j] = NULL;
    if (unrefarray (l->m[i], l->len, sizeof (char *))) {
      l->m[i] = NULL;
      if (--l->nc == 0) {
	// fixme: need global depth in parallel
	return false;
	
	free (l->m);
	free (l);
	return true; // layer has been destroyed
      }
    }
  }
  return false;
}
#endif // dimension == 3

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
  /* the current cell index and level */
  int i, j;
#if dimension == 3
  int k;
#endif
  int level;
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
#if dimension == 3
  c->p[c->n].k = p.k;
#endif
  c->n++;
}

static void cache_append (Cache * c, Point p,
			  int k, int l, int n,
			  short flags)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (Index)*c->nm);
  }
  c->p[c->n].i = p.i + k;
  c->p[c->n].j = p.j + l;
#if dimension == 3
  c->p[c->n].k = p.k + n;
#endif  
  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

@define cache_append_x(c,p,k,l,n,flags) cache_append(c,p,k,l,n,flags)
#if dimension == 2
@define cache_append_y(c,p,k,l,n,flags) cache_append(c,p,l,k,n,flags)
#else // dimension == 3
@define cache_append_y(c,p,k,l,n,flags) cache_append(c,p,n,k,l,flags)
@define cache_append_z(c,p,k,l,n,flags) cache_append(c,p,l,n,k,flags)
#endif

/* low-level memory management */
#if dimension == 2
@def allocated(k,l,n) (quadtree->L[point.level]->m[point.i+k] &&
		       quadtree->L[point.level]->m[point.i+k][point.j+l])
@
@define NEIGHBOR(k,l,n)	(quadtree->L[point.level]->m[point.i+k][point.j+l])
@def PARENT(k,l,n) (quadtree->L[point.level-1]->m[(point.i+GHOSTS)/2+k]
		  [(point.j+GHOSTS)/2+l])
@
@def allocated_child(k,l,n)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k]
   && quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k][2*point.j-GHOSTS+l])
@			   
@def CHILD(k,l,n)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k]
		 [2*point.j-GHOSTS+l])
@
#else // dimension == 3
@def allocated(a,l,n) (quadtree->L[point.level]->m[point.i+a] &&
			 quadtree->L[point.level]->m[point.i+a][point.j+l] &&
			 quadtree->L[point.level]->m[point.i+a][point.j+l]
			 [point.k+n])
@
@def NEIGHBOR(a,l,n)	(quadtree->L[point.level]->m[point.i+a][point.j+l]
			                                       [point.k+n])
@			
@def PARENT(a,l,n) (quadtree->L[point.level-1]->m[(point.i+GHOSTS)/2+a]
		  [(point.j+GHOSTS)/2+l][(point.k+GHOSTS)/2+n])
@
@def allocated_child(a,l,n)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+a]
   && quadtree->L[point.level+1]->m[2*point.i-GHOSTS+a][2*point.j-GHOSTS+l]
   && quadtree->L[point.level+1]->m[2*point.i-GHOSTS+a][2*point.j-GHOSTS+l]
			      [2*point.k-GHOSTS+n])
@
@def CHILD(a,l,n)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+a]
		    [2*point.j-GHOSTS+l][2*point.k-GHOSTS+n])
@
#endif
@define CELL(m) (*((Cell *)(m)))

/***** Multigrid macros *****/
@define depth()        (quadtree->depth)
@define aparent(k,l,n) CELL(PARENT(k,l,n))
@define child(k,l,n)   CELL(CHILD(k,l,n))

/***** Quadtree macros ****/
@define NN               (1 << point.level)
@define cell		 CELL(NEIGHBOR(0,0,0))
@define neighbor(k,l,n)	 CELL(NEIGHBOR(k,l,n))
@define neighborp(k,l,n) neighborpf(point,k+0,l+0,n+0)

static Point neighborpf (Point p, int k, int l, int n) {
#if dimension == 2
  return (Point){p.i + k, p.j + l, p.level};
#else
  return (Point){p.i + k, p.j + l, p.k + n, p.level};
#endif
}
			
/***** Data macros *****/
@define data(k,l,n)     ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
@define fine(a,k,l,n)   ((double *) (CHILD(k,l,n) + sizeof(Cell)))[a]
@define coarse(a,k,l,n) ((double *) (PARENT(k,l,n) + sizeof(Cell)))[a]

@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
#if dimension == 2
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1
  };
#else
  struct { int x, y, z; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1, 2*((point.k+GHOSTS)%2)-1
  };
#endif
  NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
  parent.j = (point.j + GHOSTS)/2;
#if dimension == 3
  parent.k = (point.k + GHOSTS)/2;
#endif
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
#define _BACK   (2*point.k - GHOSTS)
#define _FRONT  (_BACK + 1)

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
#if dimension == 2
#define _push(b,c,d,e,f)						\
  { _s++; stack[_s].l = b;						\
    stack[_s].i = c; stack[_s].j = d;					\
    stack[_s].stage = f; }
#define _pop()					  \
  { point.level = stack[_s].l;			  \
    point.i = stack[_s].i; point.j = stack[_s].j; \
    stage = stack[_s].stage; _s--; }
#else // dimension == 3
#define _push(b,c,d,e,f)						\
  { _s++; stack[_s].l = b;						\
    stack[_s].i = c; stack[_s].j = d; stack[_s].k = e;			\
    stack[_s].stage = f; }
#define _pop()								\
  { point.level = stack[_s].l;						\
    point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; \
    stage = stack[_s].stage; _s--; }
#endif

@def foreach_cell()
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 2
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    Point point = {GHOSTS,GHOSTS,GHOSTS,0};
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, GHOSTS, GHOSTS, GHOSTS, 0); /* the root cell */
    if (is_active(cell))
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0,0,0))
	continue;
      switch (stage) {
      case 0: {
	POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell()
        if (point.level < quadtree->depth) {
	  _push (point.level, point.i, point.j, point.k, 1);
          _push (point.level + 1, _LEFT, _TOP, _FRONT, 0);
        }
        break;
      }
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	      _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	      _push (point.level + 1, _LEFT,  _BOTTOM, _FRONT, 0); break;
#if dimension == 2
      case 3: _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0); break;
#else // dimension == 3
      case 3: _push (point.level, point.i, point.j, point.k, 4);
	      _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0); break;
      case 4: _push (point.level, point.i, point.j, point.k, 5);
	      _push (point.level + 1, _LEFT, _TOP, _BACK, 0); break;
      case 5: _push (point.level, point.i, point.j, point.k, 6);
	      _push (point.level + 1, _RIGHT, _TOP, _BACK, 0); break;
      case 6: _push (point.level, point.i, point.j, point.k, 7);
	      _push (point.level + 1, _LEFT, _BOTTOM, _BACK, 0); break;
      case 7: _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0); break;
#endif
      }
    }
  }
@

@def foreach_cell_post(condition)
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 2
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    Point point = {GHOSTS,GHOSTS,GHOSTS,0};
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, GHOSTS, GHOSTS, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0,0,0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
	if (point.level == quadtree->depth) {
	  _push (point.level, point.i, point.j, point.k, 8);
	}
	else {
	  _push (point.level, point.i, point.j, point.k, 1);
	  if (condition)
	    _push (point.level + 1, _LEFT, _TOP, _FRONT, 0);
	}
	break;
      }
      case 1:
	_push (point.level, point.i, point.j, point.k, 2);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0);
	break;
      case 2:
	_push (point.level, point.i, point.j, point.k, 3);
	if (condition)
	  _push (point.level + 1, _LEFT,  _BOTTOM, _FRONT, 0);
	break;
      case 3:
	_push (point.level, point.i, point.j, point.k, 4);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0);
	break;
#if dimension == 3
      case 4:
	_push (point.level, point.i, point.j, point.k, 5);
	if (condition)
	  _push (point.level + 1, _LEFT, _TOP, _BACK, 0);
	break;
      case 5:
	_push (point.level, point.i, point.j, point.k, 6);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _BACK, 0);
	break;
      case 6:
	_push (point.level, point.i, point.j, point.k, 7);
	if (condition)
	  _push (point.level + 1, _LEFT,  _BOTTOM, _BACK, 0);
	break;
      case 7:
	_push (point.level, point.i, point.j, point.k, 8);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0);
	break;	
#endif
      default: {
        POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell_post()
      }
      }
    }
  }
@

#if dimension == 2
@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 4; _k++) {
    point.i = _i + _k/2; point.j = _j + _k%2;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
}
@
#else // dimension == 3
@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  int _k = 2*point.k - GHOSTS;
  point.level++;
  for (int _l = 0; _l < 8; _l++) {
    point.i = _i + _l/4; point.j = _j + (_l/2)%2; point.k = _k + _l%2;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.k = (_k + GHOSTS)/2;
  point.level--;
}
@
#endif // dimension == 3
  
#define update_cache() { if (quadtree->dirty) update_cache_f(); }

#define is_prolongation(cell) (!is_leaf(cell) && !cell.neighbors && \
			       cell.pid >= 0)
#define is_boundary(cell) (cell.pid < 0)

@def foreach_cache(_cache,clause) {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 3
  int kg = 0; NOT_UNUSED(kg);
#endif
  OMP_PARALLEL()
#if dimension == 2
  Point point = {GHOSTS,GHOSTS,0};
#else
  Point point = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
  int _k; short _flags; NOT_UNUSED(_flags);
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
#if dimension == 3
    point.k = _cache.p[_k].k;
#endif
    point.level = _cache.p[_k].level;
    _flags = _cache.p[_k].flags;
    POINT_VARIABLES;
@
@define end_foreach_cache() } OMP_END_PARALLEL() }

@def foreach_cache_level(_cache,_l,clause) {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 3
  int kg = 0; NOT_UNUSED(kg);
#endif
  OMP_PARALLEL()
#if dimension == 2
  Point point = {GHOSTS,GHOSTS,0};
#else
  Point point = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
  point.level = _l;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
#if dimension == 3
    point.k = _cache.p[_k].k;
#endif
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
      // look in a 3x3 neighborhood for boundary cells
      // fixme: should probably be a 5x5 neighborhood
      // fixme: foreach_neighbor() would be nicer
      for (int k = -1; k <= 1; k++)
	for (int l = -1; l <= 1; l++)
#if dimension == 3
	  for (int n = -1; n <= 1; n++)
#endif
	    if (is_boundary(neighbor(k,l,n)) &&
		!(neighbor(k,l,n).flags & fboundary)) {
	      point.i += k; point.j += l;
#if dimension == 3
	      point.k += n;
#endif
	      cache_level_append (&q->boundary[level], point);
	      cell.flags |= fboundary;
	      point.i -= k; point.j -= l;
#if dimension == 3
	      point.k -= n;
#endif
	    }
    if (is_leaf (cell)) {
      if (is_local(cell)) {
	cache_append (&q->leaves, point, 0, 0, 0, 0);
	// faces
	unsigned flags = 0;
	foreach_dimension()
	  if (is_boundary(neighbor(-1)) || is_prolongation(neighbor(-1)) ||
	      is_leaf(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, 0, 0, 0, flags);
	foreach_dimension()
	  if (is_boundary(neighbor(1)) || is_prolongation(neighbor(1)) ||
	      (!is_local(neighbor(1)) && is_leaf(neighbor(1))))
	    cache_append_x (&q->faces, point, 1, 0, 0, face_x);
	// vertices: fixme 3D
	cache_append (&q->vertices, point, 0, 0, 0, 0);
	if (!is_leaf(neighbor(1,1)) || !is_local(neighbor(1,1)))
	  cache_append (&q->vertices, point, 1, 1, 0, 0);
	foreach_dimension()
	  if ((!is_leaf(neighbor(0,-1)) || !is_local(neighbor(0,-1))) &&
	      (!is_leaf(neighbor(1,0)) || !is_local(neighbor(1,0))))
	    cache_append_x (&q->vertices, point, 1, 0, 0, 0);
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
	cell.flags &= ~halo;
      }
      else if (!is_boundary(cell) || is_local(aparent(0))) { // non-local
	// faces
	unsigned flags = 0;
	foreach_dimension()
	  if (allocated(-1) &&
	      is_local(neighbor(-1)) && is_prolongation(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, 0, 0, 0, flags);
	foreach_dimension()
	  if (allocated(1) && is_local(neighbor(1)) &&
	      is_prolongation(neighbor(1)))
	    cache_append_x (&q->faces, point, 1, 0, 0, face_x);
	// vertices: fixme: this is too complicated, 3D?
	foreach_dimension()
	  for (int i = -1; i <= 1; i += 2)
	    if (allocated(i) && is_local(neighbor(i)) &&
		is_prolongation(neighbor(i))) {
	      cache_append_x (&q->vertices, point, (i + 1)/2, 0, 0, 0);
	      cache_append_x (&q->vertices, point, (i + 1)/2, 1, 0, 0);
	    }
      }
      continue;
    }
    else { // not a leaf
      bool restriction = level > 0 && (aparent(0).flags & halo);
      for (int k = -GHOSTS; k <= GHOSTS && !restriction; k++)
	for (int l = -GHOSTS; l <= GHOSTS && !restriction; l++)
#if dimension == 3
	  for (int n = -GHOSTS; n <= GHOSTS && !restriction; n++)
#endif
	    if (allocated(k,l,n) &&
		((is_local(neighbor(k,l,n)) && is_leaf(neighbor(k,l,n))) ||
		 is_remote_leaf(neighbor(k,l,n))))
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
      foreach_child()
        if (is_boundary(cell))
	  cell.flags |= halo;
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
#if dimension == 3
@define is_face_z() (_flags & face_z)
#endif

@def foreach_vertex(clause)
  update_cache();
  foreach_cache(quadtree->vertices, clause) {
    x -= Delta/2.; y -= Delta/2.;
#if dimension == 3
    z -= Delta/2.;
#endif
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
    if (is_active(cell) && is_local(cell)) {
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
#if dimension == 2
	    for (scalar s in list)
	      if (!is_constant(s))
		((double *)(L->m[i][j] + sizeof(Cell)))[s] = undefined;
#else // dimension == 3
          for (int k = 0; k < L->len; k++)
	    if (L->m[i][j][k])
	      for (scalar s in list)
  	        if (!is_constant(s))
		  ((double *)(L->m[i][j][k] + sizeof(Cell)))[s] = undefined;
#endif    
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

#if dimension == 2
static void alloc_children (Point point, int i, int j, int k)
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
      assert (!CHILD(k,l,0));
      CHILD(k,l,0) = b;
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
	alloc_children (point, k, l, 0);
      neighbor(k,l).neighbors++;
    }
}
				
static void free_children (Point point, int i, int j)
{
  point.i += i; point.j += j;
  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  assert (CHILD(0,0,0));
  mempool_free (L->pool, CHILD(0,0,0));
  for (int k = 0; k < 2; k++) {
    for (int l = 0; l < 2; l++)
      CHILD(k,l,0) = NULL;
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
#else // dimension == 3
static void alloc_children (Point point, int i, int j, int k)
{
  if (point.level == quadtree->depth)
    alloc_layer();
  point.i += i; point.j += j; point.k += k;
  
  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = mempool_alloc0 (L->pool);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      layer_add_row (L, 2*point.i - GHOSTS + k, 2*point.j - GHOSTS + l);
      for (int n = 0; n < 2; n++) {
	assert (!CHILD(k,l,n));
	CHILD(k,l,n) = b;
	b += len;
      }
    }

  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      for (int n = 0; n < 2; n++)
	child(k,l,n).pid = cell.pid;

@if TRASH
  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      for (int n = 0; n < 2; n++)
	for (scalar s in all)
	  fine(s,k,l,n) = undefined;
@endif
}

void increment_neighbors (Point point)
{
  quadtree->dirty = true;
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
      for (int n = -GHOSTS/2; n <= GHOSTS/2; n++) {
	if (neighbor(k,l,n).neighbors == 0)
	  alloc_children (point, k, l, n);
	neighbor(k,l,n).neighbors++;
      }
}

static void free_children (Point point, int i, int j, int k)
{
  point.i += i; point.j += j; point.k += k;
  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  assert (CHILD(0,0,0));
  mempool_free (L->pool, CHILD(0,0,0));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      for (int n = 0; n < 2; n++)
	CHILD(k,l,n) = NULL;
      if (layer_remove_row (L, 2*point.i - GHOSTS + k, 2*point.j - GHOSTS + l)) {
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
      for (int n = -GHOSTS/2; n <= GHOSTS/2; n++)
	if (allocated(k,l,n)) {
	  neighbor(k,l,n).neighbors--;
	  if (neighbor(k,l,n).neighbors == 0)
	    free_children (point,k,l,n);
	}
}
#endif // dimension == 3

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
#if dimension == 2
      q->L[0]->m[i][j] = realloc (q->L[0]->m[i][j], newlen);
#else
      for (int k = 0; k < len; k++)
	q->L[0]->m[i][j][k] = realloc (q->L[0]->m[i][j][k], newlen);
#endif
  /* all other levels */
  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    size_t len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << dimension)*newlen);
    for (int i = 0; i < len; i += 2)
      if (L->m[i])
	for (int j = 0; j < len; j += 2)
	  if (L->m[i][j]) {
#if dimension == 2
	    char * new = mempool_alloc (L->pool);
	    for (int k = 0; k < 2; k++)
	      for (int o = 0; o < 2; o++) {
		memcpy (new, L->m[i+k][j+o], oldlen);
		L->m[i+k][j+o] = new;
		new += newlen;
	      }
#else // dimension == 3
	    for (int k = 0; k < len; k += 2)
	      if (L->m[i][j][k]) {
		char * new = mempool_alloc (L->pool);
		for (int l = 0; l < 2; l++)
		  for (int m = 0; m < 2; m++)
		    for (int n = 0; n < 2; n++) {
		      memcpy (new, L->m[i+l][j+m][k+n], oldlen);
		      L->m[i+l][j+m][k+n] = new;
		      new += newlen;
		    }
	      }
#endif
	  }
    mempool_destroy (oldpool);
  }
}

/* Boundaries */

#define bid(cell) (- cell.pid - 1)

@define VN v.x
@define VT v.y

#define is_neighbor(...) (allocated(__VA_ARGS__) && \
			  !is_boundary(neighbor(__VA_ARGS__)) \
			  && cond(neighborp(__VA_ARGS__)))

#if _MPI
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi()  enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_coarsen (Point point, scalar s);
  
void box_boundaries (int l,
		     bool (cond1)(Point), bool (cond)(Point),
		     scalar * list)
{
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.refine != no_coarsen) {
      if (s.v.x == s) {
	if (s.face)
	  faces = vectors_add (faces, s.v);
	else
	  vectors = vectors_add (vectors, s.v);
      }
      else if (s.v.x < 0)
	scalars = list_add (scalars, s);
    }
  
  foreach_boundary (l)
    if (cond1 (point)) {
      int nc = 0, id = bid(cell);
      coord nv;
      foreach_dimension()
	nv.x = 0;
      for (scalar s in scalars)
	s[] = 0.;
      for (vector v in vectors)
	foreach_dimension()
	  v.x[] = 0.;
      // normal directions
      foreach_dimension()
	for (int i = -1; i <= 1; i += 2)
	  if (is_neighbor(i)) {
	    Point neighbor = neighborp(i);
	    for (scalar s in scalars)
	      s[] += s.boundary[id](neighbor, point, s);
	    for (vector v in vectors) {
	      scalar vn = VN;
	      v.x[] += vn.boundary[id](neighbor, point, v.x);
	    }
	    for (vector v in faces) {
	      scalar vn = VN;
	      if (v.x.normal && vn.boundary[id])
		v.x[(i + 1)/2] = vn.boundary[id](neighbor, point, v.x);
	    }
	    nc++; nv.x++;
	  }
      if (vectors || faces)
#if dimension == 2
	foreach_dimension() {
	  // tangential directions
	  if (nv.y == 0)
	    for (int i = -1; i <= 1; i += 2)
	      if (is_neighbor(i)) {
		Point neighbor = neighborp(i);
		for (vector v in vectors) {
		  scalar vt = VT;
		  v.y[] += vt.boundary[id](neighbor, point, v.y);
		}
		nv.y++;
	      }
	  // tangential face directions
	  if (faces && allocated(-1) && is_boundary(neighbor(-1))) {
	    int nf = 0;
	    for (vector v in faces)
	      v.x[] = 0.;
	    for (int j = -1; j <= 1; j += 2)
	      if (is_neighbor(0,j) || is_neighbor(-1,j)) {
		Point neighbor = neighborp(0,j);
		for (vector v in faces) {
		  scalar vt = VT;
		  v.x[] += vt.boundary[id](neighbor, point, v.x);
		}
		nf++;
	      }
	    if (nf > 0) {
	      for (vector v in faces)
		v.x[] /= nf;
	    }
	    else
	      for (vector v in faces)
		v.x[] = 0.; // fixme: this should be undefined
	  }
	}
#else // dimension == 3
      assert (false);
#endif
      // 2D diagonal directions
      if (nc == 0)
#if dimension == 3
	foreach_dimension()
#endif
	  for (int k = -1; k <= 1; k += 2)
	    if (is_boundary(neighbor(k,0))) {
	      Point n1 = neighborp(k,0);
	      int id1 = bid(neighbor(k,0));
	      for (int l = -1; l <= 1; l += 2)
		if (is_boundary(neighbor(0,l)) && is_neighbor(k,l)) {
		  Point n = neighborp(k,l), n2 = neighborp(0,l);
		  int id2 = bid(neighbor(0,l));
		  for (scalar s in scalars)
		    s[] += (s.boundary[id1](n,n1,s) + s.boundary[id2](n,n2,s) -
			    s[k,l]);
		  for (vector v in vectors) {
		    scalar vt = VT, vn = VN;
		    v.x[] += (vt.boundary[id1](n,n1,v.x) +
			      vn.boundary[id2](n,n2,v.x) -
			      v.x[k,l]);
		    v.y[] += (vn.boundary[id1](n,n1,v.y) +
			      vt.boundary[id2](n,n2,v.y) -
			      v.y[k,l]);
#if dimension == 3
		    v.z[] += (vt.boundary[id1](n,n1,v.z) +
			      vt.boundary[id2](n,n2,v.z) -
			      v.z[k,l]);
#endif
		  }
		  nv.x++; nv.y++;
#if dimension == 3
		  nv.z++;
#endif
		  nc++;
		}
	    }
      // 3D diagonal directions
#if dimension == 3
      if (nc == 0)
	for (int k = -1; k <= 1; k += 2)
	  for (int l = -1; l <= 1; l += 2)
	    if (is_boundary(neighbor(k,l,0))) {
	      Point n1 = neighborp(k,l,0);
	      int id1 = bid(neighbor(k,l,0));
	      for (int n = -1; n <= 1; n += 2)
		if (is_boundary(neighbor(k,0,n)) &&
		    is_boundary(neighbor(0,l,n)) &&
		    is_neighbor(k,l,n)) {
		  Point n0 = neighborp(k,l,n), n2 = neighborp(k,0,n);
		  Point n3 = neighborp(0,l,n);
		  int id2 = bid(neighbor(k,0,n)), id3 = bid(neighbor(0,l,n));
		  for (scalar s in scalars)
		    s[] += (s.boundary[id1](n0,n1,s) +
			    s.boundary[id2](n0,n2,s) +
			    s.boundary[id3](n0,n3,s) -
			    2.*s[k,l,n]);
		  nc++;
		}
	    }
#endif
      // averaging
      if (nc > 0) {
	for (scalar s in scalars)
	  s[] /= nc;
      }
      else
	for (scalar s in scalars)
	  s[] = undefined;
      foreach_dimension() {
	if (nv.x > 0) {
	  for (vector v in vectors)
	    v.x[] /= nv.x;
	}
	else
	  for (vector v in vectors)
	    v.x[] = undefined;
      }
    }
  free (scalars);
  free (vectors);
  free (faces);
}

#undef is_neighbor
 
@undef VN
@undef VT
@define VN _attribute[s].v.x
@define VT _attribute[s].v.y
 
static bool retrue (Point point) { return true; }

static bool retleaf (Point point) { return is_leaf(cell); }

static bool retleafhalo (Point point) {
  // leaf or prolongation or restriction halo
  return is_leaf(cell) || !cell.neighbors || (cell.flags & halo);
}

static bool retleafhalo1 (Point point) {
  // leaf or restriction halo
  return is_leaf(cell) || (cell.flags & halo);  
}

static bool rethalo (Point point) {
  // restriction halo
  return cell.flags & halo;
}
 
static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  if (l < 0)
    for (l = 0; l <= depth(); l++)
      box_boundaries (l, retrue, retleaf, list);
  else
    box_boundaries (l, retrue, retrue, list);
  enable_fpe_for_mpi();
}

static void box_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  disable_fpe_for_mpi();
  if (l == 0)
    return;
  box_boundaries (l, rethalo, retleafhalo1, list);
  enable_fpe_for_mpi();
}

static void box_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, 
					    int l, int depth)
{
  disable_fpe_for_mpi();
  if (l == depth)
    box_boundaries (l, retrue, retrue, list);
  else
    box_boundaries (l, retrue, retleafhalo, list);
  enable_fpe_for_mpi();
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
	  foreach_child() {				\
	    cell.flags &= ~(leaf|active);		\
	    cell.pid = pid;				\
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
  Point point = {0,0,0};
  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;
  int j;
  OMP(omp for schedule(static))
  for (j = 0; j < n + 2*GHOSTS; j++)
#if dimension == 3
    for (int k = 0; k < n + 2*GHOSTS; k++)
#endif
    {
      for (int i = 0; i < GHOSTS; i++)
	if (allocated(i + n,j,k))
	  for (scalar s in list1)
	    s[i,j,k] = s[i + n,j,k];
      for (int i = n + GHOSTS; i < n + 2*GHOSTS; i++)
	if (allocated(i - n,j,k))
	  for (scalar s in list1)
	    s[i,j,k] = s[i - n,j,k];
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

int refine_cell (Point point, scalar * list, int flag, Cache * refined);

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
#if dimension == 2
  for (int i = 0; i < L->len; i++) {
    for (int j = 0; j < L->len; j++)
      free (L->m[i][j]);
    free (L->m[i]);
  }
  free (L->m);
  free (L);
  /* all other levels */
  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])
	free (L->m[i]);
    mempool_destroy (L->pool);
    free (L->m);
    free (L);
  }
#else // dimension == 3
  for (int i = 0; i < L->len; i++) {
    for (int j = 0; j < L->len; j++) {
      for (int k = 0; k < L->len; k++)
	free (L->m[i][j][k]);
      free (L->m[i][j]);
    }
    free (L->m[i]);
  }
  free (L->m);
  free (L);
  /* all other levels */
  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i]) {
	for (int j = 0; j < L->len; j++)
	  if (L->m[i][j])
	    free (L->m[i][j]);
	free (L->m[i]);
      }
    mempool_destroy (L->pool);
    free (L->m);
    free (L);
  }
#endif // dimension == 3
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
#if dimension == 2
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
#else // dimension == 3
  for (int i = 0; i < L->len; i++)
    for (int j = 0; j < L->len; j++) {
      layer_add_row (L, i, j);
      for (int k = 0; k < L->len; k++)
	L->m[i][j][k] = calloc (1, sizeof(Cell) + datasize);
    }
  CELL(L->m[GHOSTS][GHOSTS][GHOSTS]).flags |= (leaf|active);
  for (int k = -GHOSTS; k <= GHOSTS; k++)
    for (int l = -GHOSTS; l <= GHOSTS; l++)
      for (int n = -GHOSTS; n <= GHOSTS; n++)
	CELL(L->m[GHOSTS+k][GHOSTS+l][GHOSTS+n]).pid =
	  (k < 0 ? -1 - left :
	   k > 0 ? -1 - right :
	   l > 0 ? -1 - top :
	   l < 0 ? -1 - bottom :
	   n > 0 ? -1 - back :
	   n < 0 ? -1 - front :
	   pid());
#endif // dimension == 3
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
#if dimension == 3
    point.k = (p.z - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS &&
	point.j >= 0 && point.j < n + 2*GHOSTS
#if dimension == 3
	&& point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_local(cell) && is_leaf(cell))
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
