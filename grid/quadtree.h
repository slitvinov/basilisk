#define GRIDNAME "Quadtree"

#include <stdio.h>
#include <assert.h>
#include "common.h"

#define I (point.i - GHOSTS)
#define J (point.j - GHOSTS)

typedef struct {
  char flags, neighbors;
  Data d;
} Cell;

typedef struct _Quadtree Point;
typedef struct _Quadtree Quadtree;

struct _Quadtree {
  int depth;       /* the maximum depth of the tree */

  Quadtree * back; /* back pointer to the "parent" quadtree */
  Cell ** m;       /* the grids at each level */
  int i, j, level; /* the current cell index and level */
};

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

/***** Quadtree macros ****/
#define _n (1 << point.level) /* fixme */
#define aparent(k,l) point.m[point.level-1][((point.i+GHOSTS)/2+k)*(_n/2+2*GHOSTS) + \
					    (point.j+GHOSTS)/2+l]
#define child(k,l)   point.m[point.level+1][(2*point.i-GHOSTS+k)*2*(_n + GHOSTS) + \
					    (2*point.j-GHOSTS+l)]
#define cell         point.m[point.level][point.i*(_n + 2*GHOSTS) + point.j]
#define parent       aparent(0,0)
#define alloc_children() { if (point.level == point.depth) alloc_layer(&point); }
#define free_children()

/***** Quadtree variables *****/
#define QUADTREE_VARIABLES \
  int    level = point.level;                                   NOT_UNUSED(level);   \
  int    childx = 2*((point.i+GHOSTS)%2)-1;                     NOT_UNUSED(childx);  \
  int    childy = 2*((point.j+GHOSTS)%2)-1;                     NOT_UNUSED(childy);  \
  double delta = 1./(1 << point.level);                         NOT_UNUSED(delta);   \

/***** Data macros *****/
#define data(k,l)  point.m[point.level][(point.i + k)*(_n + 2*GHOSTS) +	\
					(point.j + l)].d
#define fine(a,k,l) field(child(k,l).d, a, double)
#define coarse(a,k,l) field(aparent(k,l).d, a, double)

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

#if 0
void recursive (void * m, int n, int i, int j, int nl)
{
  if (n == nl) {
    /* do something */
  }
  if (n < nl) {
    m = finer_level (m, n), n *= 2;
    recursive (m, n, _LEFT,  _TOP,    nl);
    recursive (m, n, _RIGHT, _TOP,    nl);
    recursive (m, n, _LEFT,  _BOTTOM, nl);
    recursive (m, n, _RIGHT, _BOTTOM, nl);
  }
}
#endif

#define STACKSIZE 20
#define _push(b,c,d,e) _s++;						\
  { stack[_s].l = b; stack[_s].i = c; stack[_s].j = d; stack[_s].stage = e; }
#define _pop(b,c,d,e)							\
  { b = stack[_s].l; c = stack[_s].i; d = stack[_s].j; e = stack[_s].stage; _s--; }

#define foreach_boundary_cell(grid,dir)					\
  {									\
    Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid);	\
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
  }

#define foreach_cell(grid) foreach_boundary_cell(grid,0)
#define end_foreach_cell()						\
      if (point.level < point.depth) {					\
	  _push (point.level, point.i, point.j, 1);			\
          _push (point.level + 1, _LEFT, _TOP, 0);			\
      }									\
      break;								\
      }									\
      case 1: _push (point.level, point.i, point.j, 2);			\
              _push (point.level + 1, _RIGHT, _TOP,    0);	        \
      break;							\
      case 2: _push (point.level, point.i, point.j, 3);		\
              _push (point.level + 1, _LEFT,  _BOTTOM, 0); break;	\
      case 3:								\
              _push (point.level + 1, _RIGHT, _BOTTOM, 0); break;	\
      }								\
    }                                                                   \
  }

#define foreach_cell_post(grid,condition)				\
  {									\
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
  }

/* ================== derived traversals ========================= */

enum {
  active = 1 << 0,
  leaf   = 1 << 1,
  halo   = 1 << 2
};

#define foreach_leaf(grid)        foreach_cell(grid) if (cell.flags & leaf) {
#define end_foreach_leaf()        continue; } end_foreach_cell()

#define foreach       foreach_leaf
#define end_foreach   end_foreach_leaf

#define foreach_boundary(grid,dir) foreach_boundary_cell(grid,dir)
#define end_foreach_boundary()    if (cell.flags & leaf) continue; end_foreach_boundary_cell()

#define foreach_fine_to_coarse(grid) foreach_cell_post(grid,!(cell.flags & leaf))
#define end_foreach_fine_to_coarse() end_foreach_cell_post()

#define foreach_level(grid,l)      foreach_cell(grid) { \
                                      if (level == l || cell.flags & leaf) {
#define end_foreach_level()           continue; } } end_foreach_cell()

#define foreach_boundary_level(grid,dir,l)   foreach_boundary(grid,dir) {	\
                                             if (level == l || cell.flags & leaf) {
#define end_foreach_boundary_level()         continue; } } end_foreach_boundary()

void alloc_layer (Quadtree * point)
{
  Quadtree * q = point->back;
  q->depth++; point->depth++;
  q->m = &(q->m[-1]);
  q->m = realloc(q->m, sizeof (Cell *)*(q->depth + 2)); 
  q->m = &(q->m[1]);
  point->m = q->m;
  q->m[q->depth] = calloc (_size(q->depth), sizeof (Cell));
}

/* fixme: this needs to be merged with refine_wavelet() */
int refine_quadtree (Quadtree * quadtree)
{
  int nf = 0;
  foreach_leaf (quadtree) {
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
      }
    nf++;
  } end_foreach_leaf();
  return nf;
}

void * init_grid (int n)
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
  q->depth = 0;
  q->m = malloc(sizeof (Cell *)*2);
  q->m[0] = NULL; q->m = &(q->m[1]); /* make sure we don't try to access level -1 */
  q->m[0] = calloc (_size(0), sizeof (Cell));
  /* initialise the root cell */
  q->m[0][2 + 2*GHOSTS].flags |= (leaf | active);
  q->m[0][2 + 2*GHOSTS].neighbors = 1; // only itself as neighbor
  while (depth--)
    refine_quadtree(q);
  return (void *) q;
}

void free_grid (void * m)
{
  Quadtree * q = m;
  for (int l = 0; l <= q->depth; l++)
    free (q->m[l]);
  q->m = &(q->m[-1]);
  free(q->m);
  free(q);
}
