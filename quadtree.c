#define GRID "Quadtree"

#include <stdio.h>
#include <assert.h>
#include "utils.h"

#define I (_q.i - GHOSTS)
#define J (_q.j - GHOSTS)

typedef struct {
  char flags, neighbors;
  Data d;
} Cell;

typedef struct {
  Cell ** m;       /* the grids at each level */
  int depth;       /* the maximum depth of the tree */
  int i, j, level; /* the current cell index and level */
} Quadtree;

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

int mlevel (int n)
{
  int r = 0;
  while (n > 1) { n /= 2; r++; }
  return r;
}

/***** Quadtree macros ****/
#define _n (1 << _q.level)
#define aparent(k,l) _q.m[_q.level-1][((_q.i+GHOSTS)/2+k)*(_n/2+2*GHOSTS) + (_q.j+GHOSTS)/2+l]
#define child(k,l)   _q.m[_q.level+1][(2*_q.i-GHOSTS+k)*2*(_n + GHOSTS) + (2*_q.j-GHOSTS+l)]
#define cell         _q.m[_q.level][_q.i*(_n + 2*GHOSTS) + _q.j]
#define parent       aparent(0,0)

/***** Quadtree variables *****/
#define QUADTREE_VARIABLES \
  int  level = _q.level;                                   NOT_UNUSED(level);   \
  int  childx = 2*((_q.i+GHOSTS)%2)-1;                     NOT_UNUSED(childx);  \
  int  childy = 2*((_q.j+GHOSTS)%2)-1;                     NOT_UNUSED(childy);

/***** Data macros *****/
#define data(k,l)  _q.m[_q.level][(_q.i + k)*(_n + 2*GHOSTS) + (_q.j + l)].d
#define fine(a,k,l) field(child(k,l).d, a, double)
#define coarse(a,k,l) field(aparent(k,l).d, a, double)
#define foreach_fine(a,b) for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)

/* ===============================================================
 *                    Quadtree traversal
 * traverse_recursive() below is for reference only. The macro
 * foreach_leaf() is a stack-based implementation of
 * traverse_recursive(). It is about 12% slower than the recursive
 * version and 60% slower than simple array traversal.
 *
 * This article was useful:
 * http://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
 *
 * =============================================================== */

#define _BOTTOM (2*_q.j - GHOSTS)
#define _TOP    (_BOTTOM + 1)
#define _LEFT   (2*_q.i - GHOSTS)
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

void traverse_recursive (void * m, int n)
{
  m = (void *) (((char *)m) +  sizeof(Cell)*(_totalsize(mlevel(n)) - _size(0)));
  recursive (m, 1, GHOSTS, GHOSTS, n); /* level 0, central cell */
}
#endif

#define STACKSIZE 20
#define push(b,c,d,e) _s++;						\
  { stack[_s].l = b; stack[_s].i = c; stack[_s].j = d; stack[_s].stage = e; }
#define pop(b,c,d,e)							\
  { b = stack[_s].l; c = stack[_s].i; d = stack[_s].j; e = stack[_s].stage; _s--; }

#define foreach_boundary_cell(m,n,dir)					\
  {									\
    Quadtree _q = *((Quadtree *)m);					\
    int _d = dir; NOT_UNUSED(_d);					\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      pop (_q.level, _q.i, _q.j, stage);				\
      switch (stage) {							\
      case 0: {								\
        QUADTREE_VARIABLES;						\
	VARIABLES;							\
	/* do something */
#define end_foreach_boundary_cell()					\
        if (_q.level < _q.depth) {					\
	  push (_q.level, _q.i, _q.j, 1);				\
	  int k = _d > left ? _LEFT : _RIGHT - _d;			\
	  int l = _d < top  ? _TOP  : _TOP + 2 - _d;			\
          push (_q.level + 1, k, l, 0);					\
        }								\
	break;								\
      }								        \
      case 1: {								\
	  int k = _d > left ? _RIGHT : _RIGHT - _d;			\
	  int l = _d < top  ? _BOTTOM  : _TOP + 2 - _d;			\
          push (_q.level + 1, k, l, 0);			                \
          break;                                                        \
        }								\
      }									\
    }                                                                   \
  }

#define foreach_cell(m,n) foreach_boundary_cell(m,n,0)
#define end_foreach_cell()						\
      if (_q.level < _q.depth) {					\
	  push (_q.level, _q.i, _q.j, 1);					\
          push (_q.level + 1, _LEFT, _TOP, 0);				       \
        }								       \
	break;								       \
      }								               \
      case 1: push (_q.level, _q.i, _q.j, 2); push (_q.level + 1, _RIGHT, _TOP,    0); break; \
      case 2: push (_q.level, _q.i, _q.j, 3); push (_q.level + 1, _LEFT,  _BOTTOM, 0); break; \
      case 3:                                 push (_q.level + 1, _RIGHT, _BOTTOM, 0); break; \
      }								               \
    }                                                                          \
  }

#define foreach_cell_post(m,n,condition)				\
  {									\
    Quadtree _q = *((Quadtree *)m);					\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      pop (_q.level, _q.i, _q.j, stage);				\
      switch (stage) {							\
      case 0: {								\
        QUADTREE_VARIABLES;						\
	VARIABLES;							\
	if (condition) {						\
	  if (_q.level == _q.depth)	{				\
	    push (_q.level, _q.i, _q.j, 4);				\
	  }								\
	  else {							\
	    push (_q.level, _q.i, _q.j, 1);				\
	    push (_q.level + 1, _LEFT, _TOP, 0);			\
	  }								\
	}								\
	break;								       \
      }								               \
      case 1: push (_q.level, _q.i, _q.j, 2); push (_q.level + 1, _RIGHT, _TOP,    0); break; \
      case 2: push (_q.level, _q.i, _q.j, 3); push (_q.level + 1, _LEFT,  _BOTTOM, 0); break; \
      case 3: push (_q.level, _q.i, _q.j, 4); push (_q.level + 1, _RIGHT, _BOTTOM, 0); break; \
      case 4: {								       \
        QUADTREE_VARIABLES;						\
        VARIABLES;							\
	/* do something */
#define end_foreach_cell_post()						\
      }									       \
      }								               \
    }                                                                          \
  }

/* ================== derived traversals ========================= */

enum {
  leaf     = 1 << 0,
  inactive = 1 << 1,
  halo     = 1 << 2
};

#define foreach_leaf(m,n)         foreach_cell(m,n) if (cell.flags & leaf) {
#define end_foreach_leaf()        continue; } end_foreach_cell()

#define foreach       foreach_leaf
#define end_foreach   end_foreach_leaf

#define foreach_boundary(m,n,dir) foreach_boundary_cell(m,n,dir)
#define end_foreach_boundary()    if (cell.flags & leaf) continue; end_foreach_boundary_cell()

#define foreach_fine_to_coarse(m,n) foreach_cell_post(m,n,!(cell.flags & leaf))
#define end_foreach_fine_to_coarse() end_foreach_cell_post()

void * init_grid (int n)
{
  Quadtree * q = malloc(sizeof (Quadtree));
  q->depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      free(q);
      exit (1);
    }
    n /= 2;
    q->depth++;
  }
  q->m = malloc(sizeof (Cell *)*(q->depth + 2));
  q->m[0] = malloc(sizeof(Cell)); q->m = &(q->m[1]); /* add a dummy grid at level -1 */
  for (int l = 0; l <= q->depth; l++)
    q->m[l] = calloc (_size(l), sizeof (Cell));
  Quadtree _q = *q;
  for (_q.level = 0; _q.level < _q.depth; _q.level++)
    for (_q.i = GHOSTS; _q.i < (1 << _q.level) + GHOSTS; _q.i++)
      for (_q.j = GHOSTS; _q.j < (1 << _q.level) + GHOSTS; _q.j++)
	cell.neighbors = (2*GHOSTS + 1)*(2*GHOSTS + 1);
  for (_q.i = GHOSTS; _q.i < (1 << _q.level) + GHOSTS; _q.i++)
    for (_q.j = GHOSTS; _q.j < (1 << _q.level) + GHOSTS; _q.j++) {
      cell.flags |= leaf;
      cell.neighbors = (2*GHOSTS + 1)*(2*GHOSTS + 1);
    }
  return (void *) q;
}

void free_grid (void * m)
{
  Quadtree * q = m;
  for (int l = 0; l <= q->depth; l++)
    free (q->m[l]);
  q->m--;
  free(q->m[0]);
  free(q->m);
  free(q);
}
