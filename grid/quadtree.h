#define GRIDNAME "Quadtree"

#include <stdio.h>
#include <assert.h>

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

typedef struct {
  char flags, neighbors;
} Cell;

typedef struct _Quadtree Point;
typedef struct _Quadtree Quadtree;

struct _Quadtree {
  int depth;       /* the maximum depth of the tree */

  Quadtree * back; /* back pointer to the "parent" quadtree */
  char ** m;       /* the grids at each level */
  int i, j, level; /* the current cell index and level */
};

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))

/***** Multigrid macros *****/
#define depth()      (((Quadtree *)grid)->depth)
#define aparent(k,l) CELL(point.m, point.level-1, ((point.i+GHOSTS)/2+k)*(_n/2+2*GHOSTS) + \
			  (point.j+GHOSTS)/2+l)
#define child(k,l)   CELL(point.m, point.level+1, (2*point.i-GHOSTS+k)*2*(_n + GHOSTS) + \
			  (2*point.j-GHOSTS+l))

/***** Quadtree macros ****/
#define _n (1 << point.level) /* fixme */
#define cell               CELL(point.m, point.level, point.i*(_n + 2*GHOSTS) + point.j)
#define parent             aparent(0,0)
#define alloc_children()   { if (point.level == point.depth) alloc_layer(&point); }
#define free_children()

/***** Quadtree variables *****/
#define QUADTREE_VARIABLES \
  int    level = point.level;                                   NOT_UNUSED(level);   \
  int    childx = 2*((point.i+GHOSTS)%2)-1;                     NOT_UNUSED(childx);  \
  int    childy = 2*((point.j+GHOSTS)%2)-1;                     NOT_UNUSED(childy);

/***** Data macros *****/
#define data(k,l)  ((double *) &point.m[point.level][((point.i + k)*(_n + 2*GHOSTS) + \
			       (point.j + l))*(sizeof(Cell) + datasize) + sizeof(Cell)])
#define field(cell) ((double *)(((char *) &cell) + sizeof(Cell)))
#define fine(a,k,l) field(child(k,l))[a]
#define coarse(a,k,l) field(aparent(k,l))[a]

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
    QUADTREE_VARIABLES;
    VARIABLES;
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
#define _push(b,c,d,e)						\
  { _s++; stack[_s].l = b; stack[_s].i = c; stack[_s].j = d; stack[_s].stage = e; }
#define _pop(b,c,d,e)							\
  { b = stack[_s].l; c = stack[_s].i; d = stack[_s].j; e = stack[_s].stage; _s--; }

#define foreach_boundary_cell(dir)					\
  OMP_PARALLEL()							\
    Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid); \
    int _d = dir; NOT_UNUSED(_d);					\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
        QUADTREE_VARIABLES;						\
	VARIABLES;							\
	/* do something */
#define end_foreach_boundary_cell()					\
      if (point.level < point.depth) {					\
	  _push (point.level, point.i, point.j, 1);			\
	  int k = _d > left ? _LEFT : _RIGHT - _d;			\
	  int l = _d < top  ? _TOP  : _TOP + 2 - _d;			\
          _push (point.level + 1, k, l, 0);				\
        }								\
	break;								\
      }								        \
      case 1: {								\
	  int k = _d > left ? _RIGHT : _RIGHT - _d;			\
	  int l = _d < top  ? _BOTTOM  : _TOP + 2 - _d;			\
          _push (point.level + 1, k, l, 0);				\
          break;                                                        \
        }								\
      }									\
    }                                                                   \
  OMP_END_PARALLEL()

#define foreach_cell() foreach_boundary_cell(0)
#define end_foreach_cell()						\
      if (point.level < point.depth) {					\
	  _push (point.level, point.i, point.j, 1);			\
          _push (point.level + 1, _LEFT, _TOP, 0);			\
      }									\
      break;								\
      }									\
      case 1: _push (point.level, point.i, point.j, 2);			\
              _push (point.level + 1, _RIGHT, _TOP,    0); break;	\
      case 2: _push (point.level, point.i, point.j, 3);		        \
              _push (point.level + 1, _LEFT,  _BOTTOM, 0); break;	\
      case 3: _push (point.level + 1, _RIGHT, _BOTTOM, 0); break;	\
      }								        \
    }                                                                   \
  OMP_END_PARALLEL()

#define foreach_cell_post(condition)					\
  OMP_PARALLEL()							\
    Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid);	\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
        QUADTREE_VARIABLES;						\
	VARIABLES;							\
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
        QUADTREE_VARIABLES;						\
        VARIABLES;							\
	/* do something */
#define end_foreach_cell_post()						\
      }									\
      }								        \
    }                                                                   \
  OMP_END_PARALLEL()

/* ================== derived traversals ========================= */

enum {
  active = 1 << 0,
  leaf   = 1 << 1,
  halo   = 1 << 2
};

#define foreach_leaf()            foreach_cell() if (cell.flags & leaf) {
#define end_foreach_leaf()        continue; } end_foreach_cell()

#define foreach()     foreach_leaf()
#define end_foreach   end_foreach_leaf

#define foreach_boundary(dir)     foreach_boundary_cell(dir)
#define end_foreach_boundary()    if (cell.flags & leaf) continue; end_foreach_boundary_cell()

#define foreach_fine_to_coarse()     foreach_cell_post(!(cell.flags & leaf))
#define end_foreach_fine_to_coarse() end_foreach_cell_post()

#define foreach_level(l)            foreach_cell() { \
                                      if (level == l || cell.flags & leaf) {
#define end_foreach_level()           continue; } } end_foreach_cell()

#define foreach_boundary_level(dir,l)      foreach_boundary(dir) {	\
                                             if (level == l || cell.flags & leaf) {
#define end_foreach_boundary_level()         continue; } } end_foreach_boundary()

void alloc_layer (Quadtree * p)
{
  Quadtree * q = p->back;
  q->depth++; p->depth++;
  q->m = &(q->m[-1]);
  q->m = realloc(q->m, sizeof (char *)*(q->depth + 2)); 
  q->m = &(q->m[1]);
  p->m = q->m;
  q->m[q->depth] = calloc (_size(q->depth), sizeof (Cell) + datasize);
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

Point refine_cell (Point point, scalar start, scalar end);

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
  Quadtree * q = malloc(sizeof (Quadtree));
  q->depth = 0; q->i = q->j = GHOSTS; q->level = 0.;
  q->m = malloc(sizeof (char *)*2);
  q->m[0] = NULL; q->m = &(q->m[1]); /* make sure we don't try to access level -1 */
  /* initialise the root cell */
  q->m[0] = calloc (_size(0), sizeof (Cell) + datasize);
  CELL(q->m, 0, 2 + 2*GHOSTS).flags |= (leaf | active);
  CELL(q->m, 0, 2 + 2*GHOSTS).neighbors = 1; // only itself as neighbor
  grid = q;
  while (depth--)
    foreach_leaf ()
      point = refine_cell (point, 0, nvar - 1);
}

void free_grid (void)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++)
    free (q->m[l]);
  q->m = &(q->m[-1]);
  free(q->m);
  free(q);
}

// The functions below should be independent from the details of the implementation

Point refine_cell (Point point, scalar start, scalar end)
{
  QUADTREE_VARIABLES;
  alloc_children();
  cell.flags &= ~leaf;
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      assert(!(child(k,l).flags & active));
      child(k,l).flags |= (active | leaf);
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors++;
      /* bilinear interpolation from coarser level */
      for (scalar v = start; v <= end; v++)
	fine(v,k,l) = 
	  (9.*val(v,0,0) + 3.*(val(v,2*k-1,0) + val(v,0,2*l-1)) + val(v,2*k-1,2*l-1))/16.;
    }
  return point;
}

Point locate (double xp, double yp)
{
  foreach_cell () {
    double delta = DELTA;
    double x = (point.i - GHOSTS + 0.5)*delta - 0.5;
    double y = (point.j - GHOSTS + 0.5)*delta - 0.5;
    delta /= 2.;
    if (xp < x - delta || xp > x + delta || yp < y - delta || yp > y + delta)
      continue;
    if (cell.flags & leaf)
      return point;
  }
  Point point = {-1, NULL, NULL, -1, -1, -1}; // not found
  return point;
}
