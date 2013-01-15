#include <stdlib.h>
#include <assert.h>

typedef struct _Data Data;

/* multigrid quadtree implementation */

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

int level (int n)
{
  int r = 0;
  while (n > 1) { n /= 2; r++; }
  return r;
}

void * mgrid (int r, size_t s)
{
  return malloc (s*totalsize (r));
}

#define fine(k,l) mf[(2*i-1+k)*(2*n + 2) + (2*j-1+l)]

Data * coarser_level (Data * m, int n)
{
  return (Data *) (((char *)m) + sizeof(Data)*(n+2)*(n+2));
}

Data * finer_level (Data * m, int n)
{
  return (Data *) (((char *)m) - sizeof(Data)*(2*n+2)*(2*n+2));
}

#define foreach(m) for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++)

#define foreach_fine_to_coarse(m,n)		   \
  Data * mf = m; m = coarser_level (m, n), n /= 2; \
  for (int l = level(n); \
       n > 0; \
       mf = m, m = coarser_level (m, n), n /= 2, l--) \
    for (int i = 1; i <= n; i++) \
      for (int j = 1; j <= n; j++)

/* ===============================================================
 *                    Quadtree traversal
 * traverse_recursive() below is for reference only. The macro
 * foreach_leaf() is a stack-based implementation of
 * traverse_recursive(). It is about 12% slower than the recursive
 * version and 60% slower than simple array traversal.
 *
 * This article was useful:
 * http://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
 */

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
  m = (Data *) (((char *)m) +  sizeof(Data)*(totalsize(level(n)) - 9));
  recursive (m, 1, 1, 1, n); /* level 0 */
}

#define STACKSIZE 20
#define push(b,c,d,e) top++; \
  stack[top].l = b; stack[top].i = c; stack[top].j = d; stack[top].stage = e;
#define pop(b,c,d,e) \
  b = stack[top].l; c = stack[top].i; d = stack[top].j; e = stack[top].stage; \
  top--;

#define foreach_leaf() \
  { \
    int depth = level (n); \
    assert (depth < STACKSIZE); \
    struct { int l, i, j, stage; } stack[STACKSIZE]; int top = -1; /* the stack */ \
    Data * ml[STACKSIZE]; /* pointers on each level */		\
    for (int i = depth; i >= 0; i--) { \
      ml[i] = m; m = coarser_level (m, 1 << i); \
    } \
    push (0, 1, 1, 0); /* the root */ \
    while (top >= 0) { \
      int l, i, j, stage; \
      pop (l, i, j, stage); \
      switch (stage) { \
      case 0: \
	if (l == depth) { \
	  m = ml[l]; int n = 1 << l;		\
	  /* do something */
#define end_foreach_leaf() \
	} \
	if (l < depth) { push (l, i, j, 1); push (l + 1, 2*i-1, 2*j, 0); } \
	break; \
      case 1: push (l, i, j, 2); push (l + 1, 2*i,   2*j,   0); break; \
      case 2: push (l, i, j, 3); push (l + 1, 2*i-1, 2*j-1, 0); break; \
      case 3:                    push (l + 1, 2*i,   2*j-1, 0); break; \
      } \
    } \
  }
