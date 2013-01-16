#ifdef CARTESIAN
  #define GRID "Multigrid quadtree (Cartesian)"
#else
  #define GRID "Multigrid quadtree"
#endif

#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>

typedef struct _Data Data;

int size (int l)
{
  int n = (1 << l) + 2;
  return n*n;
}

int totalsize (int l)
{
  l++;
  return ((1 << (2*l)) - 1)/3 + 4*((1 << l) - 1) + 4*l;
}

int mlevel (int n)
{
  int r = 0;
  while (n > 1) { n /= 2; r++; }
  return r;
}

#define field(data,offset,type) (*((type *)(((char *) &(data)) + (offset))))
typedef size_t var;
#define var(a) offsetof(Data,a)
#define cell(k,l) _m[(i+k)*(_n + 2) + (j+l)]
#define stencil(a,k,l) field(cell(k,l), a, double)
#define val(a) stencil(a,0,0)

Data * coarser_level (Data * m, int n)
{
  return (Data *) (((char *)m) + sizeof(Data)*(n+2)*(n+2));
}

Data * finer_level (Data * m, int n)
{
  return (Data *) (((char *)m) - sizeof(Data)*(2*n+2)*(2*n+2));
}

#define foreach_level(m,n) {				\
  Data * _m = m; int _n = n;				\
  for (int i = 1; i <= _n; i++)				\
    for (int j = 1; j <= _n; j++)
#define end_foreach_level() }

#define foreach_fine_to_coarse(m,n) {				\
  Data * _m = m, * _mf = _m;					\
  _m = coarser_level (_m, n); int _n = (n)/2;			\
  for (int level = mlevel(_n);					\
       _n > 0;							\
       _mf = _m, _m = coarser_level (_m, _n), _n /= 2, level--)	\
    for (int i = 1; i <= _n; i++)				\
      for (int j = 1; j <= _n; j++)
#define end_foreach_fine_to_coarse() }

#define fine(a,k,l) field(_mf[(2*i-1+k)*(2*_n + 2) + (2*j-1+l)], a, double)
#define coarse(a,k,l) stencil(a,k,l)

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

void recursive (Data * m, int n, int i, int j, int nl)
{
  if (n == nl) {
    /* do something */
  }
  if (n < nl) {
    m = finer_level (m, n), n *= 2;
    recursive (m, n, 2*i-1, 2*j,   nl);
    recursive (m, n, 2*i,   2*j,   nl);
    recursive (m, n, 2*i-1, 2*j-1, nl);
    recursive (m, n, 2*i,   2*j-1, nl);
  }
}

void traverse_recursive (Data * m, int n)
{
  m = (Data *) (((char *)m) +  sizeof(Data)*(totalsize(mlevel(n)) - 9));
  recursive (m, 1, 1, 1, n); /* level 0 */
}

#define STACKSIZE 20
#define push(b,c,d,e) top++;						\
  stack[top].l = b; stack[top].i = c; stack[top].j = d; stack[top].stage = e;
#define pop(b,c,d,e)							\
  b = stack[top].l; c = stack[top].i; d = stack[top].j; e = stack[top].stage; \
  top--;

#define NOT_UNUSED(a) (a = a)

enum {
  leaf   = 1 << 0
};

#define foreach_cell(m,n)						\
  {									           \
    int depth = mlevel (n);						           \
    assert (depth < STACKSIZE);						           \
    struct { int l, i, j, stage; } stack[STACKSIZE]; int top = -1; /* the stack */ \
    Data * _m = m, * ml[STACKSIZE]; /* pointers on each level */	           \
    for (int i = depth; i >= 0; i--) {					           \
      ml[i] = _m; _m = coarser_level (_m, 1 << i);			           \
    }									           \
    push (0, 1, 1, 0); /* the root */					           \
    while (top >= 0) {							\
      int level, i, j, stage;						\
      pop (level, i, j, stage);						\
      switch (stage) {							\
      case 0: {								\
        Data * _m = ml[level]; NOT_UNUSED(_m);				\
        int _n = 1 << level; NOT_UNUSED(_n);				\
	/* do something */
#define end_foreach_cell()						\
        if (!(cell(0,0).flags & leaf)) {				\
          push (level, i, j, 1);					\
          push (level + 1, 2*i-1, 2*j, 0);				       \
        }								       \
	break;								       \
      }								               \
      case 1: push (level, i, j, 2); push (level + 1, 2*i,   2*j,   0); break; \
      case 2: push (level, i, j, 3); push (level + 1, 2*i-1, 2*j-1, 0); break; \
      case 3:                        push (level + 1, 2*i,   2*j-1, 0); break; \
      }								               \
    }                                                                          \
  }

/* ================== derived traversals ========================= */

#define foreach_leaf(m,n)   foreach_cell(m,n) if (cell(0,0).flags & leaf) {
#define end_foreach_leaf()  continue; } end_foreach_cell()

void * mgrid (int r, size_t s)
{
  void * m = calloc (totalsize(r), s);
  foreach_level (m, 1 << r)
    cell(0,0).flags |= leaf;
  end_foreach_level();
  return m;
}

#ifdef CARTESIAN /* uses finest grid traversal only (i.e. a purely Cartesian grid) */
  #define foreach foreach_level
  #define end_foreach end_foreach_level
#else
#define foreach foreach_leaf
#define end_foreach end_foreach_leaf
#endif
