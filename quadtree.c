#define GRID "Quadtree"

#include <stdio.h>
#include <assert.h>

typedef struct {
  int flags;
  Data d;
} Cell;

#define GHOSTS 1        // number of ghost layers
#define I (i - GHOSTS)
#define J (j - GHOSTS)

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

size_t _totalsize (int l)
{
  size_t s = 0;
  while (l >= 0)
    s += _size(l--);
  return s;
}

int mlevel (int n)
{
  int r = 0;
  while (n > 1) { n /= 2; r++; }
  return r;
}

#define CELL(p)   ((Cell *)p)
#define cell      CELL(_m)[i*(_n + 2*GHOSTS) + j]
#define data(k,l) CELL(_m)[(i + k)*(_n + 2*GHOSTS) + (j + l)].d

void * coarser_level (void * m, int n)
{
  return (void *) (((char *)m) + sizeof(Cell)*(n + 2*GHOSTS)*(n + 2*GHOSTS));
}

void * finer_level (void * m, int n)
{
  return (void *) (((char *)m) - sizeof(Cell)*4*(n + GHOSTS)*(n + GHOSTS));
}

#define finecell(k,l) CELL(_mf)[(2*i-GHOSTS+k)*2*(_n + GHOSTS) + (2*j-GHOSTS+l)]
#define fine(a,k,l) field(finecell(k,l).d, a, double)
#define foreach_fine(a,b) for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)
#define coarse stencil


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

#define _BOTTOM (2*j - GHOSTS)
#define _TOP    (_BOTTOM + 1)
#define _LEFT   (2*i - GHOSTS)
#define _RIGHT  (_LEFT + 1)

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

#define STACKSIZE 20
#define push(b,c,d,e) _s++;						\
  stack[_s].l = b; stack[_s].i = c; stack[_s].j = d; stack[_s].stage = e;
#define pop(b,c,d,e)							\
  b = stack[_s].l; c = stack[_s].i; d = stack[_s].j; e = stack[_s].stage; _s--;

#define foreach_boundary_cell(m,n,dir)					\
  {									           \
    int depth = mlevel (n), _d = dir; NOT_UNUSED(_d);				   \
    assert (depth < STACKSIZE);						           \
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    void * _m = m, * _ml[STACKSIZE]; /* pointers on each level */	           \
    for (int i = depth; i >= 0; i--) {					           \
      _ml[i] = _m; _m = coarser_level (_m, 1 << i);			           \
    }									           \
    push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int level, i, j, stage;						\
      pop (level, i, j, stage);						\
      switch (stage) {							\
      case 0: {								\
        void * _m = _ml[level]; NOT_UNUSED(_m);				\
        int _n = 1 << level; NOT_UNUSED(_n);				\
	VARIABLES;							\
	/* do something */
#define end_foreach_boundary_cell()					\
        if (level < depth) {						\
          push (level, i, j, 1);					\
	  int k = _d > left ? _LEFT : _RIGHT - _d;			\
	  int l = _d < top  ? _TOP  : _TOP + 2 - _d;			\
          push (level + 1, k, l, 0);					\
        }								\
	break;								\
      }								        \
      case 1: {								\
	  int k = _d > left ? _RIGHT : _RIGHT - _d;			\
	  int l = _d < top  ? _BOTTOM  : _TOP + 2 - _d;			\
          push (level + 1, k, l, 0);			                \
          break;                                                        \
        }								\
      }									\
    }                                                                   \
  }

#define foreach_cell(m,n) foreach_boundary_cell(m,n,0)
#define end_foreach_cell()						\
        if (level < depth) {						\
          push (level, i, j, 1);					\
          push (level + 1, _LEFT, _TOP, 0);				       \
        }								       \
	break;								       \
      }								               \
      case 1: push (level, i, j, 2); push (level + 1, _RIGHT, _TOP,    0); break; \
      case 2: push (level, i, j, 3); push (level + 1, _LEFT,  _BOTTOM, 0); break; \
      case 3:                        push (level + 1, _RIGHT, _BOTTOM, 0); break; \
      }								               \
    }                                                                          \
  }

#define foreach_cell_post(m,n,condition)				\
  {									\
    int depth = mlevel (n);						\
    assert (depth < STACKSIZE);						           \
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */  \
    void * _m = m, * _ml[STACKSIZE]; /* pointers on each level */	           \
    for (int i = depth; i >= 0; i--) {					           \
      _ml[i] = _m; _m = coarser_level (_m, 1 << i);			           \
    }									           \
    push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int level, i, j, stage;						\
      pop (level, i, j, stage);						\
      switch (stage) {							\
      case 0: {								\
        void * _m = _ml[level]; NOT_UNUSED(_m);				\
        int _n = 1 << level;    NOT_UNUSED(_n);				\
	VARIABLES;							\
	if (level < depth && (condition)) {				\
          push (level, i, j, 1);					\
          push (level + 1, _LEFT, _TOP, 0);				       \
        }								       \
	break;								       \
      }								               \
      case 1: push (level, i, j, 2); push (level + 1, _RIGHT, _TOP,    0); break; \
      case 2: push (level, i, j, 3); push (level + 1, _LEFT,  _BOTTOM, 0); break; \
      case 3: push (level, i, j, 4); push (level + 1, _RIGHT, _BOTTOM, 0); break; \
      case 4: {								       \
        void * _m = _ml[level];    NOT_UNUSED(_m);				\
        void * _mf = _ml[level+1]; NOT_UNUSED(_mf);				\
        int _n = 1 << level;       NOT_UNUSED(_n);				\
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
  inactive = 1 << 1
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
  int r = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    r++;
  }
  void * _m = calloc (_totalsize(r), sizeof (Cell));
  int _n = 1 << r;
  for (int i = GHOSTS; i < _n + GHOSTS; i++)
    for (int j = GHOSTS; j < _n + GHOSTS; j++)
      cell.flags |= leaf;
  return _m;
}

#define free_grid free
