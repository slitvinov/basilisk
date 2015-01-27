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

@if _MPI
@ define is_local(cell)  ((cell).pid == pid())
@else
@ define is_local(cell)  true
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

// Layer

typedef struct {
  char *** m; // the 2D array of data
  int * nr;   // the number of allocated rows for each column
  int nc;     // the number of allocated columns
  int len;    // the (1D) size of the array
} Layer;

static size_t _size (size_t l)
{
  return (1 << l) + 2*GHOSTS;
}

static Layer * new_layer (int depth)
{
  Layer * l = malloc (sizeof (Layer));
  l->len = _size (depth);
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
  CacheLevel * active;   /* active cells indices for each level */
  CacheLevel * prolongation; /* halo prolongation indices for each level */
  CacheLevel * restriction;  /* halo restriction indices for each level */

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

/* low-level memory management */
@define allocated(k,l) (quadtree->L[point.level]->m[point.i+k] && quadtree->L[point.level]->m[point.i+k][point.j+l])
@define NEIGHBOR(k,l)	(quadtree->L[point.level]->m[point.i+k][point.j+l])
@define PARENT(k,l) (quadtree->L[point.level-1]->m[(point.i+GHOSTS)/2+k][(point.j+GHOSTS)/2+l])
@define CHILD(k,l)  (quadtree->L[point.level+1]->m[2*point.i-GHOSTS+k][2*point.j-GHOSTS+l])
@define CELL(m) (*((Cell *)(m)))

/***** Multigrid macros *****/
@define depth()      (quadtree->depth)
@define aparent(k,l) CELL(PARENT(k,l))
@define child(k,l)   CELL(CHILD(k,l))

/***** Quadtree macros ****/
@define NN              (1 << point.level)
@define cell		CELL(NEIGHBOR(0,0))
@define neighbor(k,l)	CELL(NEIGHBOR(k,l))

/***** Data macros *****/
@define data(k,l)     ((double *) (NEIGHBOR(k,l) + sizeof(Cell)))
@define fine(a,k,l)   ((double *) (CHILD(k,l) + sizeof(Cell)))[a]
@define coarse(a,k,l) ((double *) (PARENT(k,l) + sizeof(Cell)))[a]

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
    Point point = {GHOSTS,GHOSTS,0};
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
        if (point.level < quadtree->depth) {
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

#define update_cache() { if (quadtree->dirty) update_cache_f(); }

static void update_cache_f (void)
{
  Quadtree * q = grid;

  /* empty caches */
  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n = q->restriction[l].n = 0;

  foreach_cell() {
    if (is_active(cell) && is_local(cell)) { // always true in serial
      // active cells
      cache_level_append (&q->active[level], point);
      if (is_leaf (cell)) {
	cache_append (&q->leaves, point, 0, 0, 0);
	// faces
	unsigned flags = 0;
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
	bool restriction =
	  level > 0 &&
	  is_local(aparent(0,0)) &&
	  (aparent(0,0).flags & halo);
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
    else if (!is_active(cell))
      continue;
  }
  
  q->dirty = false;
}

@def foreach_cache(_cache,clause) {
  update_cache();
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

@define foreach(clause) foreach_cache(quadtree->leaves, clause)
@define end_foreach()   end_foreach_cache()

@def foreach_face_generic(clause) 
  foreach_cache(quadtree->faces, clause) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() (_flags & face_x)
@define is_face_y() (_flags & face_y)

@def foreach_vertex(clause) 
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
}

static void alloc_children (Point point, int i, int j)
{
  if (point.level == quadtree->depth)
    alloc_layer();
  point.i += i; point.j += j;

  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = calloc (4, len);
  for (int k = 0; k < 2; k++) {
    layer_add_row (L, 2*point.i - GHOSTS + k);
    for (int l = 0; l < 2; l++) {
      assert (!CHILD(k,l));
      CHILD(k,l) = b;
      b += len;
    }
  }

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
}

static void free_children (Point point, int i, int j)
{
  point.i += i; point.j += j;
  /* low-level memory management */
  Layer * L = quadtree->L[point.level + 1];
  assert (CHILD(0,0));
  free (CHILD(0,0));
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
    char *** m = q->L[l]->m;
    size_t len = q->L[l]->len;
    for (int i = 0; i < len; i += 2)
      if (m[i])
	for (int j = 0; j < len; j += 2)
	  if (m[i][j]) {
	    char * new = malloc (4*newlen), * old = m[i][j];
	    for (int k = 0; k < 2; k++)
	      for (int o = 0; o < 2; o++) {
		memcpy (new, m[i+k][j+o], oldlen);
		m[i+k][j+o] = new;
		new += newlen;
	      }
	    free (old);
	  }
  }
}

@def foreach_halo(name,_l) {
  update_cache();
  CacheLevel _cache = quadtree->name[_l];
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

void refine_cell (Point point, scalar * list, int flag, int * nactive);

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
      if (L->m[i]) {
	if (i%2 == 0)
	  for (int j = 0; j < L->len; j += 2)
	    free (L->m[i][j]);
	free (L->m[i]);
      }
    free (L->m);
    free (L->nr);
    free (L);
  }
  q->L = &(q->L[-1]);
  free (q->L);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->restriction);
  free (q);
  grid = NULL;
}

void init_grid (int n)
{
  // check 64 bits structure alignment
  assert (sizeof(Cell)%8 == 0);
  
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
@if _MPI
  for (int k = -GHOSTS; k <= GHOSTS; k++)
    for (int l = -GHOSTS; l <= GHOSTS; l++)
      CELL(L->m[GHOSTS+k][GHOSTS+l]).pid = -1;
  CELL(L->m[GHOSTS][GHOSTS]).pid = pid();
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
      refine_cell (point, NULL, 0, NULL);
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
