#define GRIDNAME "Quadtree"
#define dimension 2

#define DYNAMIC 1 // use dynamic data allocation
#define TWO_ONE 1 // enforce 2:1 refinement ratio
#define GHOSTS  2

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

typedef struct {
  unsigned flags;
  unsigned neighbors; // number of refined neighbors in a 3x3 neighborhood
@if _MPI
  int pid; // fixme: alignment ?
@endif
} Cell;

enum {
  active  = 1 << 0,
  leaf    = 1 << 1,
  refined = 1 << 2,
  halo    = 1 << 3,
  border  = 1 << 4,
  user    = 5,

  face_x = 1 << 0,
  face_y = 1 << 1
};

#define _CORNER 4
#define is_active(cell)  ((cell).flags & active)
#define is_leaf(cell)    ((cell).flags & leaf)
#define is_refined(cell) (is_active(cell) && !is_leaf(cell))
#define is_corner(cell)  (stage == _CORNER)
#define is_coarse()      (!is_leaf(cell))
#define is_border(cell)  ((cell).flags & border)

// Caches

typedef struct {
  int i, j;
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
  int i, j;
  short level, flags;
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

// Quadtree

typedef struct _Point Quadtree;
struct _Point {
  int depth;        /* the maximum depth of the tree */

  Quadtree * back;  /* back pointer to the "parent" quadtree */
#if DYNAMIC
  char *** m;       /* the grids at each level */
#else
  char ** m;        /* the grids at each level */
#endif
  int i, j, level;  /* the current cell index and level */

  Cache        leaves;   /* leaf indices */
  Cache        faces;    /* face indices */
  Cache        vertices; /* vertex indices */
  CacheLevel * active;   /* active cells indices for each level */
  CacheLevel * prolongation; /* halo prolongation indices for each level */
  CacheLevel * restriction;  /* halo restriction indices for each level */

  bool dirty;       /* whether caches should be updated */
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

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#if DYNAMIC
  @define CELL(m,level,i) (*((Cell *)m[level][i]))
  @define allocated(i,j) (point.m[point.level][_index(i,j)])
#else
  @define CELL(m,level,i) (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))
  @define allocated(i,j) true
#endif

/***** Multigrid macros *****/
@define depth()      (((Quadtree *)grid)->depth)
@define _index(k,l)  ((point.i + k)*(NN + 2*GHOSTS) + point.j + l)
@def _parentindex(k,l) (((point.i+GHOSTS)/2+k)*(NN/2+2*GHOSTS) +
			(point.j+GHOSTS)/2+l) @
@def _childindex(k,l) ((2*point.i-GHOSTS+k)*2*(NN + GHOSTS) +
		       (2*point.j-GHOSTS+l)) @
@define aparent(k,l) CELL(point.m, point.level-1, _parentindex(k,l))
@define child(k,l)   CELL(point.m, point.level+1, _childindex(k,l))

/***** Quadtree macros ****/
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

@def foreach_cell()
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Quadtree point = *((Quadtree *)grid); point.back = grid;
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
        if (point.level < point.depth) {
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
    Quadtree point = *((Quadtree *)grid); point.back = grid;
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
	if (point.level == point.depth) {
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

#define corners()							\
      if (_corners) {							\
        if (_d < top) {							\
  	  if (point.j == GHOSTS)					\
	    _push (point.level, point.i, point.j - 1, _CORNER);		\
	  if (point.j == GHOSTS + NN - 1)			        \
	    _push (point.level, point.i, point.j + 1, _CORNER);		\
	} else {							\
	  if (point.i == GHOSTS)					\
	    _push (point.level, point.i - 1, point.j, _CORNER);		\
	  if (point.i == GHOSTS + NN - 1)			        \
	    _push (point.level, point.i + 1, point.j, _CORNER);		\
        }								\
      }

@def foreach_boundary_cell(dir,corners)
  { _OMPSTART /* for face reduction */
    int ig = _ig[dir], jg = _jg[dir];	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Quadtree point = *((Quadtree *)grid); point.back = grid;
    int _d = dir; NOT_UNUSED(_d);
    /* traverse corners only for top and bottom */
    bool _corners = (corners);
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop (point.level, point.i, point.j, stage);
      if (!allocated (0,0))
	continue;
      switch (stage) {
      case 0: case _CORNER: {
	POINT_VARIABLES;
	/* do something */
@
@def end_foreach_boundary_cell()
        }
	if (is_corner (cell)) continue;
        /* children */
        if (point.level < point.depth) {
	  _push (point.level, point.i, point.j, 1);
	  int k = _d > left ? _LEFT : _RIGHT - _d;
	  int l = _d < top  ? _TOP  : _TOP + 2 - _d;
	  _push (point.level + 1, k, l, 0);
	} else corners();
	break;
      case 1: {
  	  int k = _d > left ? _RIGHT : _RIGHT - _d;
	  int l = _d < top  ? _BOTTOM  : _TOP + 2 - _d;
	  _push (point.level + 1, k, l, 0);
	  corners();
	  break;
        }
      }
    }  _OMPEND
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
@define end_foreach_child_direction() end_foreach_child()

#define update_cache() { if (((Quadtree *)grid)->dirty) update_cache_f(); }

static void update_cache_f (void)
{
  Quadtree * q = grid;

  /* empty caches */
  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n = q->restriction[l].n = 0;

  foreach_cell()
    if (is_active(cell)) { // always true in serial
      // active cells
      cache_level_append (&q->active[level], point);
      if (is_leaf (cell)) {
	cache_append (&q->leaves, point, 0, 0, 0);
	// faces
	short flags = 0;
	foreach_dimension()
	  if (!is_refined(neighbor(-1,0)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, 0, 0, flags);
	if (!is_active(neighbor(1,0)))
	  cache_append (&q->faces, point, 1, 0, face_x);
	if (!is_active(neighbor(0,1)))
	  cache_append (&q->faces, point, 0, 1, face_y);
	// vertices
	cache_append (&q->vertices, point, 0, 0, 0);
	if (!is_leaf(neighbor(1,1)))
	  cache_append (&q->vertices, point, 1, 1, 0);
	if (!is_leaf(neighbor(0,-1)) && !is_leaf(neighbor(1,0)))
	  cache_append (&q->vertices, point, 1, 0, 0);
	if (!is_leaf(neighbor(-1,0)) && !is_leaf(neighbor(0,1)))
	  cache_append (&q->vertices, point, 0, 1, 0);
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
	cell.flags &= ~halo;
	continue;
      }
      else { // not a leaf
	bool restriction = level > 0 ? (aparent(0,0).flags & halo) : false;
	for (int k = -GHOSTS; k <= GHOSTS && !restriction; k++)
	  for (int l = -GHOSTS; l <= GHOSTS && !restriction; l++)
	    if (allocated(k,l) && is_leaf(neighbor(k,l)))
	      restriction = true;
	if (restriction) {
	  // halo restriction
	  cell.flags |= halo;
	  cache_level_append (&q->restriction[level], point);
	}
	else
	  cell.flags &= ~halo;
      }
    }

  q->dirty = false;
}

@def foreach_cache(_cache,clause) {
  update_cache();
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Quadtree point = *((Quadtree *)grid); point.back = grid;
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
  Quadtree point = *((Quadtree *)grid); point.back = grid;
  point.level = _l;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
    POINT_VARIABLES;
@
@define end_foreach_cache_level() } OMP_END_PARALLEL() }

@define foreach(clause) foreach_cache(((Quadtree *)grid)->leaves, clause)
@define end_foreach()   end_foreach_cache()

@def foreach_face_generic(clause) 
  foreach_cache(((Quadtree *)grid)->faces, clause) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() (_flags & face_x)
@define is_face_y() (_flags & face_y)

@def foreach_vertex(clause) 
  foreach_cache(((Quadtree *)grid)->vertices, clause) {
    x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_cache()

@def foreach_fine_to_coarse(clause)     {
  update_cache();
  for (int _l = depth() - 1; _l >= 0; _l--) {
    CacheLevel _active = ((Quadtree *)grid)->active[_l];
    foreach_cache_level (_active,_l,clause)
      if (!is_leaf (cell)) {
@
@define end_foreach_fine_to_coarse() } end_foreach_cache_level(); } }

@def foreach_level(l) {
  update_cache();
  CacheLevel _active = ((Quadtree *)grid)->active[l];
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
  return m;
}
#endif // !DYNAMIC

@def cache_level_resize(name)
q->name = realloc (q->name, (q->depth + 1)*sizeof (CacheLevel));
cache_level_init (&q->name[q->depth]);
@

void alloc_layer (Quadtree * p)
{
  Quadtree * q = p->back;
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
  cache_level_resize (active);
  cache_level_resize (prolongation);
  cache_level_resize (restriction);
}

static void alloc_children (Quadtree * q, int k, int l)
{
  if (q->level == q->depth) alloc_layer(q);
  Point point = *((Point *)q);
  point.i += k; point.j += l;

#if DYNAMIC
  char ** m = ((char ***)q->m)[q->level+1];
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      if (!m[_childindex(k,l)])
	m[_childindex(k,l)] = calloc (1, sizeof(Cell) + datasize);
#endif // !DYNAMIC

@if _MPI
  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      child(k,l).pid = cell.pid;
  if (is_border(cell)) {
    bool neighbors = false;
    for (int k = -GHOSTS/2; k <= GHOSTS/2 && !neighbors; k++)
      for (int l = -GHOSTS/2; l <= GHOSTS/2 && !neighbors; l++)
	neighbors = (neighbor(k,l).pid != cell.pid);
    if (neighbors)
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++)
	  child(k,l).flags |= border;
  }
@endif
  
@if TRASH
  // foreach child
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      for (scalar s in all)
	fine(s,k,l) = undefined;
@endif
}

void increment_neighbors (Quadtree * p)
{
  p->back->dirty = true;
  Point point = *((Point *)p);
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++) {
      if (neighbor(k,l).neighbors == 0)
	alloc_children (&point, k, l);
      neighbor(k,l).neighbors++;
    }
  *((Point *)p) = point;
}

static void free_children (Point point, int k, int l)
{
#if DYNAMIC
  Quadtree * q = &point;
  point.i += k; point.j += l;
  char ** m = ((char ***)q->m)[q->level+1];
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      free (m[_childindex(k,l)]);
      m[_childindex(k,l)] = NULL;
    }
#endif
}

void decrement_neighbors (Point point)
{
  ((Quadtree *)&point)->back->dirty = true;
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
  Quadtree * q = grid;
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

@def foreach_halo(name,_l) {
  update_cache();
  CacheLevel _cache = ((Quadtree *)grid)->name[_l];
  foreach_cache_level (_cache, _l,)
@
@define end_foreach_halo() end_foreach_cache_level(); }

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * lleft, * lright;
  list_split (list, d, &lleft, &lright);

  if (l < 0)
    foreach_boundary_cell (d, true)
      if (is_leaf (cell)) {
	for (scalar s in lright) {
	  s[ghost] = s.boundary[d] (point, s);
	  point.i -= ig; point.j -= jg;
	  double vb = s.boundary[d] (point, s);
	  point.i += ig; point.j += jg;
	  s[2*ig,2*jg] = vb;
	}
	for (scalar s in lleft)
	  s[] = s.boundary[d] (point, s);
	corners(); /* we need this otherwise we'd skip corners */
	continue;
      }
  else
    foreach_boundary_cell (d, true) {
      if (level == l) {
	if (allocated(ig,jg))
	  for (scalar s in lright) {
	    s[ghost] = s.boundary[d] (point, s);
	    point.i -= ig; point.j -= jg;
	    double vb = s.boundary[d] (point, s);
	    point.i += ig; point.j += jg;
	    s[2*ig,2*jg] = vb;
	  }
	for (scalar s in lleft)
	  s[] = s.boundary[d] (point, s);
	corners(); /* we need this otherwise we'd skip corners */
	continue;
      }
      else if (is_leaf (cell))
	continue;
    }

  free (lleft);
  free (lright);
}

static void box_boundary_halo_prolongation_normal (const Boundary * b,
						   scalar * list, 
						   int l, int depth)
{
  // see test/boundary_halo.c
  int d = ((BoxBoundary *)b)->d, in, jn;
  if (d % 2)
    in = jn = 0;
  else {
    in = _ig[d]; jn = _jg[d];
  }

  foreach_boundary_cell (d, false) {
    if (level == l) {
      if (l == depth ||          // target level
	  is_leaf(cell) ||       // leaves
	  (cell.flags & halo) || // restriction halo
	  is_corner(cell)) {     // corners
	// leaf or halo restriction
	for (scalar s in list)
	  s[in,jn] = s.boundary[d] (point, s);
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d) {
	  if (allocated(in,jn))
	    for (scalar s in list)
	      s[in,jn] = s.boundary[d] (point, s);
	}
      continue;
    }
  }
}

static void box_boundary_halo_prolongation_tangent (const Boundary * b,
						    scalar * list, 
						    int l, int depth)
{
  // see test/boundary_halo.c
  int d = ((BoxBoundary *)b)->d;

  foreach_boundary_cell (d, true) {
    if (level == l) {
      if (l == depth ||          // target level
	  is_leaf(cell) ||       // leaves
	  (cell.flags & halo) || // restriction halo
	  is_corner(cell)) {     // corners
	// leaf or halo restriction
	for (scalar s in list)
	  s[ghost] = s.boundary[d] (point, s);
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d) {
	  if (allocated(ig,jg))
	    for (scalar s in list)
	      s[ghost] = s.boundary[d] (point, s);
	}
      continue;
    }
  }  
}

static void box_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, 
					    int l, int depth)
{
  // see test/boundary_halo.c
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL, * tangent = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[d]) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  if (s.normal)
	    normal = list_add (normal, s);
	}
	else
	  tangent = list_add (tangent, s);
      }	
      else
	centered = list_add (centered, s);
    }

  foreach_boundary_cell (d, d > left) {
    if (level == l) {
      if ((l == depth ||          // target level
	   is_leaf(cell) ||       // leaves
	   (cell.flags & halo) || // restriction halo
	   is_corner(cell)) &&   // corners
	  allocated(ig,jg)) {
	// leaf or halo restriction
	for (scalar s in centered) {
	  s[ghost] = s.boundary[d] (point, s);
	  point.i -= ig; point.j -= jg;
	  double vb = s.boundary[d] (point, s);
	  point.i += ig; point.j += jg;
	  s[2*ig,2*jg] = vb;
	}
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d) {
	  if (allocated(ig,jg))
	    for (scalar s in centered)
	      s[ghost] = s.boundary[d] (point, s);
	  if (allocated(2*ig,2*jg))
	    for (scalar s in centered) {
	      point.i -= ig; point.j -= jg;
	      double vb = s.boundary[d] (point, s);
	      point.i += ig; point.j += jg;
	      s[2*ig,2*jg] = vb;
	    }
	}
      continue;
    }
  }  
  free (centered);

  box_boundary_halo_prolongation_normal (b, normal, l, depth);
  free (normal);
  box_boundary_halo_prolongation_tangent (b, tangent, l, depth);
  free (tangent);
}

Point refine_cell (Point point, scalar * list, int flag, int * nactive);

static void free_cache (CacheLevel * c)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++)
    if (c[l].n > 0)
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
  for (int l = 0; l <= q->depth; l++) {
#if DYNAMIC
    for (int i = 0; i < _size(l); i++)
      free (q->m[l][i]);
#endif
    free (q->m[l]);
  }
  q->m = &(q->m[-1]);
  free(q->m);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->restriction);
  free(q);
  grid = NULL;
}

void init_grid (int n)
{
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
  CELL(q->m, 0, 2*GHOSTS*(GHOSTS + 1)).flags |= (leaf|active);
@if _MPI
  for (int k = -GHOSTS; k <= GHOSTS; k++)
    for (int l = -GHOSTS; l <= GHOSTS; l++)
      CELL(q->m, 0, ((GHOSTS + k)*(1 + 2*GHOSTS) + GHOSTS + l)).pid = -1;
  CELL(q->m, 0, 2*GHOSTS*(GHOSTS + 1)).pid = pid();
@endif
  cache_init (&q->leaves);
  cache_init (&q->faces);
  cache_init (&q->vertices);
  q->active = calloc (1, sizeof (CacheLevel));
  q->prolongation = calloc (1, sizeof (CacheLevel));
  q->restriction = calloc (1, sizeof (CacheLevel));
  q->dirty = true;
  grid = q;
  N = 1 << depth;
  while (depth--)
    foreach_leaf()
      point = refine_cell (point, NULL, 0, NULL);
@if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
@endif
  update_cache();
  trash (all);
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level = b->restriction = box_boundary_level;
    b->halo_restriction  = no_halo_restriction;
    b->halo_prolongation = box_boundary_halo_prolongation;
    add_boundary (b);
  }
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

#include "quadtree-common.h"

@if _MPI
#include "quadtree-mpi.h"
@endif
