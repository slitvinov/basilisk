#define GRIDNAME "Binary tree"
#define dimension 1

#define DYNAMIC 1 // use dynamic data allocation
#define TWO_ONE 1 // enforce 2:1 refinement ratio
#define GHOSTS  2

#define I     (point.i - GHOSTS)
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
  remote_leaf = 1 << 4,
  user    = 5
};

@define is_active(cell)  ((cell).flags & active)
@define is_leaf(cell)    ((cell).flags & leaf)
@define is_corner(cell)  false
@define is_coarse()      ((cell).neighbors > 0)
@define is_local(cell)   true
@define is_remote_leaf(cell) false

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
  int i;
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

// Bitree

typedef struct {
  int depth;        /* the maximum depth of the tree */

#if DYNAMIC
  char *** m;       /* the grids at each level */
#else
  char ** m;        /* the grids at each level */
#endif

  Cache        leaves;  /* leaf indices */
  Cache        faces;   /* face indices */
  Cache        refined;  /* refined cells */
  CacheLevel * active;  /* active cells indices for each level */
  CacheLevel * prolongation; /* halo prolongation indices for each level */
  CacheLevel * restriction;  /* halo restriction indices for each level */
  CacheLevel   coarsened; /* coarsened cells */
  
  bool dirty;       /* whether caches should be updated */
} Bitree;

#define bitree ((Bitree *)grid)
#define quadtree bitree

struct _Point {
  int i, level;  /* the current cell index and level */
};
static Point last_point;

static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (IndexLevel)*c->nm);
  }
  c->p[c->n].i = p.i;
  c->n++;
}

static void cache_append (Cache * c, Point p, int k, short flags)
{
  if (c->n >= c->nm) {
    c->nm += 100;
    c->p = realloc (c->p, sizeof (Index)*c->nm);
  }
  c->p[c->n].i = p.i + k;
  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n;
}

#if DYNAMIC
@define CELL(m,level,i) (*((Cell *)m[level][i]))
  @define allocated(i,j,k) (bitree->m[point.level][_index(i)])
  @define allocated_child(i,j,k) (bitree->m[point.level+1][_childindex(i)])
#else
  @define CELL(m,level,i) (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))
  @define allocated(i,j,k) true
  @define allocated_child(i,j,k) true
#endif

/***** Multigrid macros *****/
@define depth()      (((Bitree *)grid)->depth)
@define _index(k)    (point.i + k)
@define _parentindex(k) ((point.i + GHOSTS)/2 + k)
@define _childindex(k) (2*point.i - GHOSTS + k)
@define aparent(k,l,d) CELL(bitree->m, point.level-1, _parentindex(k))
@define child(k,l,d)   CELL(bitree->m, point.level+1, _childindex(k))

/***** Bitree macros ****/
@define cell		CELL(bitree->m, point.level, _index(0))
@define neighbor(k,l,d)	CELL(bitree->m, point.level, _index(k))

/***** Data macros *****/
#if DYNAMIC
  @def data(k,l,d)
    ((double *) (bitree->m[point.level][_index(k)] + sizeof(Cell))) @
  @def fine(a,k,l,d)
    ((double *) (bitree->m[point.level+1][_childindex(k)] + sizeof(Cell)))[a] @
  @def coarse(a,k,l,d)
    ((double *) (bitree->m[point.level-1][_parentindex(k)] + sizeof(Cell)))[a] @
#else // !DYNAMIC
  @def data(k,l,d)
    ((double *) &bitree->m[point.level]
     [_index(k)*(sizeof(Cell) + datasize) + sizeof(Cell)]) @
  @def fine(a,k,l,d)
    ((double *) &bitree->m[point.level+1]
     [_childindex(k)*(sizeof(Cell) + datasize) + sizeof(Cell)])[a] @
  @def coarse(a,k,l,d)
    ((double *) &bitree->m[point.level-1]
     [_parentindex(k)*(sizeof(Cell) + datasize) + sizeof(Cell)])[a] @
#endif // !DYNAMIC

@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
  struct { int x; } child = { 2*((point.i+GHOSTS)%2)-1 }; NOT_UNUSED(child);
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
  if (point.level == bitree->depth) {
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
    Point point = {GHOSTS,0};
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
        if (point.level < bitree->depth) {
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
	if (point.level == point.depth) {
	  _push (point.level, point.i, 4);
	}
	else {
	  _push (point.level, point.i, 1);
	  if (condition)
	    _push (point.level + 1, _LEFT, 0);
	}
	break;
      }
      case 1:
	_push (point.level, point.i, 4);
	if (condition)
	  _push (point.level + 1, _RIGHT, 0);
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

// ghost cell coordinates for each direction
static int _ig[] = {1,-1};

@define corners()
@def foreach_boundary_cell(dir,corners)
  if (dir <= left) { _OMPSTART /* for face reduction */
    int ig = _ig[dir]; NOT_UNUSED(ig);
    Point point = {GHOSTS,0};
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
        if (point.level < bitree->depth)
	  _push (point.level + 1, _d == left ? _LEFT : _RIGHT, 0);
	break;
      }
    }  _OMPEND
  }
@

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

@def foreach_child_direction(d) {
  int _i = 2*point.i - GHOSTS;
  point.level++;
  {
    point.i = _i + (d == right);
    POINT_VARIABLES;
@
@def end_foreach_child_direction()
  }
  point.i = (_i + GHOSTS)/2;
  point.level--;
}
@

#define update_cache() { if (((Bitree *)grid)->dirty) update_cache_f(); }

#define is_prolongation(cell) (!is_leaf(cell) && !cell.neighbors && \
			       cell.pid >= 0)
#define is_boundary(cell) (cell.pid < 0)

static void update_cache_f (void)
{
  Bitree * q = grid;

  /* empty caches */
  q->leaves.n = q->faces.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n = q->restriction[l].n = 0;

  foreach_cell() {
    if (is_active(cell)) { // always true in serial
      // active cells
      cache_level_append (&q->active[level], point);
      if (is_leaf (cell)) {
	cache_append (&q->leaves, point, 0, 0);
	// faces
	cache_append (&q->faces, point, 0, 0);
	if (!is_active(neighbor(1,0)))
	  cache_append (&q->faces, point, 1, 0);
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
	cell.flags &= ~halo;
	continue;
      }
      else { // not a leaf
	bool restriction = level > 0 ? (aparent(0,0).flags & halo) : false;
	for (int k = -GHOSTS; k <= GHOSTS && !restriction; k++)
	  if (allocated(k,0) && is_leaf(neighbor(k,0)))
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
@if _MPI
    // !is_active(cell)
    else if (is_leaf(cell)) {
      // non-local halo prolongation
      if (cell.neighbors > 0)
	cache_level_append (&q->prolongation[level], point);
      continue;
    }
@endif // _MPI
  }

  q->dirty = false;
}

@def foreach_cache(_cache,clause) {
  update_cache();
  int ig = 0; NOT_UNUSED(ig);
  OMP_PARALLEL()
  Point point = {GHOSTS,0};
  int _k; short _flags; NOT_UNUSED(_flags);
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.level = _cache.p[_k].level;
    _flags = _cache.p[_k].flags;
    POINT_VARIABLES;
@
@define end_foreach_cache() } OMP_END_PARALLEL() }

@def foreach_cache_level(_cache,_l,clause) {
  int ig = 0; NOT_UNUSED(ig);
  OMP_PARALLEL()
  Point point = {GHOSTS,0};
  point.level = _l;
  int _k;
  OMP(omp for schedule(static) clause)
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    POINT_VARIABLES;
@
@define end_foreach_cache_level() } OMP_END_PARALLEL() }

@define foreach(clause) foreach_cache(((Bitree *)grid)->leaves, clause)
@define end_foreach()   end_foreach_cache()

@def foreach_face_generic(clause) 
  foreach_cache(((Bitree *)grid)->faces, clause) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() (true)

@def foreach_vertex(clause)
  foreach_cache(((Bitree *)grid)->faces, clause) {
    x -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_cache()

@def foreach_fine_to_coarse(clause)     {
  update_cache();
  for (int _l = depth() - 1; _l >= 0; _l--) {
    CacheLevel _active = ((Bitree *)grid)->active[_l];
    foreach_cache_level (_active,_l,clause)
      if (!is_leaf (cell)) {
@
@define end_foreach_fine_to_coarse() } end_foreach_cache_level(); } }

@def foreach_level(l) {
  update_cache();
  CacheLevel _active = ((Bitree *)grid)->active[l];
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

@define foreach_leaf()            foreach_cell() if (is_leaf (cell)) {
@define end_foreach_leaf()        continue; } end_foreach_cell()

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
	  if (!is_constant(s))
	    ((double *)(q->m[l][i] + sizeof(Cell)))[s] = undefined;
#else // !DYNAMIC
  int cellsize = sizeof(Cell) + datasize;;
  for (int l = 0; l <= q->depth; l++) {
    char * data = q->m[l] + sizeof(Cell);
    for (int i = 0; i < _size(l); i++, data += cellsize)
      for (scalar s in list)
	if (!is_constant(s))
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

void alloc_layer (void)
{
  Bitree * q = bitree;
  q->depth++;
  q->m = &(q->m[-1]);
#if DYNAMIC
  q->m = realloc(q->m, sizeof (char **)*(q->depth + 2)); 
  q->m = &(q->m[1]);
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

void alloc_children (Point point, int k, int l)
{
  if (level == bitree->depth)
    alloc_layer();
  point.i += k;

#if DYNAMIC
  char ** m = ((char ***)bitree->m)[level+1];
  for (int k = 0; k < 2; k++)
    if (!m[_childindex(k)])
      m[_childindex(k)] = calloc (1, sizeof(Cell) + datasize);
#endif // !DYNAMIC

  // foreach child
  for (int k = 0; k < 2; k++)
    child(k).pid = cell.pid;
  
@if TRASH
  // foreach child
  for (int k = 0; k < 2; k++)
    for (scalar s in all)
      fine(s,k,l) = undefined;
@endif
}

void increment_neighbors (Point point)
{
  bitree->dirty = true;
  int l = 0;
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++) {
    if (neighbor(k,l).neighbors == 0)
      alloc_children (point, k, l);
    neighbor(k,l).neighbors++;
  }
}

static void free_children (Point point, int k, int l)
{
#if DYNAMIC
  point.i += k;
  char ** m = ((char ***)bitree->m)[level+1];
  for (int k = 0; k < 2; k++) {
    free (m[_childindex(k)]);
    m[_childindex(k)] = NULL;
  }
#endif
}

void decrement_neighbors (Point point)
{
  bitree->dirty = true;
  int l = 0;
  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++) {
    neighbor(k,l).neighbors--;
    if (neighbor(k,l).neighbors == 0)
      free_children(point,k,l);
  }
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

@def foreach_halo(name,_l) {
  update_cache();
  CacheLevel _cache = ((Bitree *)grid)->name[_l];
  foreach_cache_level (_cache, _l,)
@
@define end_foreach_halo() end_foreach_cache_level(); }

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  if (d > left)
    return;
  scalar * lleft, * lright;
  list_split (list, d, &lleft, &lright);

  if (l < 0) {
    foreach_boundary_cell (d, true)
      if (is_leaf (cell)) {
	Point neighbor = {point.i + ig}, n2 = {point.i + 2*ig};
	for (scalar s in lright) {
	  s[ig] = s.boundary[d] (point, neighbor, s);
	  point.i -= ig;
	  double vb = s.boundary[d] (point, n2, s);
	  point.i += ig;
	  s[2*ig] = vb;
	}
	for (scalar s in lleft)
	  s[] = s.boundary[d] (point, neighbor, s);
	corners(); /* we need this otherwise we'd skip corners */
	continue;
      }
  }
  else
    foreach_boundary_cell (d, true) {
      if (level == l) {
	Point neighbor = {point.i + ig}, n2 = {point.i + 2*ig};
	for (scalar s in lright) {
	  s[ig] = s.boundary[d] (point, neighbor, s);
	  point.i -= ig;
	  double vb = s.boundary[d] (point, n2, s);
	  point.i += ig;
	  s[2*ig] = vb;
	}
	for (scalar s in lleft)
	  s[] = s.boundary[d] (point, neighbor, s);
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
  int d = ((BoxBoundary *)b)->d, in;
  if (d > left)
    return;
  if (d % 2)
    in = 0;
  else {
    in = _ig[d];
  }

  foreach_boundary_cell (d, false) {
    if (level == l) {
      if (l == depth ||          // target level
	  is_leaf(cell) ||       // leaves
	  (cell.flags & halo) || // restriction halo
	  is_corner(cell)) {     // corners
	// leaf or halo restriction
	Point neighbor = {point.i + in, 0};
	for (scalar s in list)
	  s[in,jn] = s.boundary[d] (point, neighbor, s);
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d)
	  if (allocated(in,jn)) {
	    Point neighbor = {point.i + in, 0};
	    for (scalar s in list)
	      s[in,jn] = s.boundary[d] (point, neighbor, s);
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
  if (d > left)
    return;

  foreach_boundary_cell (d, true) {
    if (level == l) {
      if (l == depth ||          // target level
	  is_leaf(cell) ||       // leaves
	  (cell.flags & halo) || // restriction halo
	  is_corner(cell)) {     // corners
	// leaf or halo restriction
	Point neighbor = {point.i + ig};
	for (scalar s in list)
	  s[ig] = s.boundary[d] (point, neighbor, s);
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d)
	  if (allocated(ig,jg)) {
	    Point neighbor = {point.i + ig};
	    for (scalar s in list)
	      s[ig] = s.boundary[d] (point, neighbor, s);
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
  if (d > left)
    return;
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
      if (l == depth ||          // target level
	  is_leaf(cell) ||       // leaves
	  (cell.flags & halo) || // restriction halo
	  is_corner(cell)) {     // corners
	// leaf or halo restriction
	Point neighbor = {point.i + ig}, n2 = {point.i + 2*ig};
	for (scalar s in centered) {
	  s[ig] = s.boundary[d] (point, neighbor, s);
	  point.i -= ig;
	  double vb = s.boundary[d] (point, n2, s);
	  point.i += ig;
	  s[2*ig] = vb;
	}
	corners();
      }
      continue;
    }
    else if (is_leaf(cell)) {
      if (level == l - 1 && cell.neighbors > 0)
	// halo prolongation
	foreach_child_direction(d) {
	  Point neighbor = {point.i + ig}, n2 = {point.i + 2*ig};
	  if (allocated(ig,jg))
	    for (scalar s in centered)
	      s[ig] = s.boundary[d] (point, neighbor, s);
	  if (allocated(2*ig,2*jg))
	    for (scalar s in centered) {
	      point.i -= ig;
	      double vb = s.boundary[d] (point, n2, s);
	      point.i += ig;
	      s[2*ig] = vb;
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
/* Periodic boundaries */

static double periodic_bc (Point, Point, scalar);
  
static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  Point point = {0,0};
  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;
  for (int i = 0; i < GHOSTS; i++)
    if (allocated(i + n,0))
      for (scalar s in list1)
	s[i,0] = s[i + n,0];
  for (int i = n + GHOSTS; i < n + 2*GHOSTS; i++)
    if (allocated(i - n,0))
      for (scalar s in list1)
	s[i,0] = s[i - n,0];
  
  free (list1);
}

static void periodic_boundary_halo_prolongation_x (const Boundary * b,
						   scalar * list, 
						   int l, int depth)
{
  periodic_boundary_level_x (b, list, l);
}
 
void refine_cell (Point point, scalar * list, int flag, Cache * refined);

static void free_cache (CacheLevel * c)
{
  Bitree * q = grid;
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
  Bitree * q = grid;
  free (q->leaves.p);
  free (q->faces.p);
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
  q = calloc (1, sizeof (Bitree));
  q->depth = 0;
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
  CELL(q->m, 0, GHOSTS).flags |= (leaf|active);
  for (int k = -GHOSTS; k <= GHOSTS; k++)
    CELL(q->m, 0, GHOSTS+k).pid = -1;
  CELL(q->m, 0, GHOSTS).pid = pid();
  cache_init (&q->leaves);
  cache_init (&q->faces);
  q->active = calloc (1, sizeof (CacheLevel));
  q->prolongation = calloc (1, sizeof (CacheLevel));
  q->restriction = calloc (1, sizeof (CacheLevel));
  q->dirty = true;
  grid = q;
  N = 1 << depth;
  while (depth--)
    foreach_leaf()
      refine_cell (point, NULL, 0, NULL);
  update_cache();
  trash (all);
@if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
@endif
  // box boundaries
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = calloc (1, sizeof (BoxBoundary));
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level = b->restriction = box_boundary_level;
    b->halo_restriction  = no_halo_restriction;
    b->halo_prolongation = box_boundary_halo_prolongation;
    add_boundary (b);
  }
  // periodic boundaries
  Boundary * b = calloc (1, sizeof (Boundary));
  b->level = b->restriction = periodic_boundary_level_x;
  b->halo_prolongation = periodic_boundary_halo_prolongation_x;
  add_boundary (b);
  init_events();
}

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    double a = (p.x - X0)/L0*n;
    if (a >= 0.5 - GHOSTS && a < n + GHOSTS - 0.5) {
      point.i = a + GHOSTS;
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

void bitree_methods()
{
  quadtree_methods();
}
