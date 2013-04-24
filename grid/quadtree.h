#define GRIDNAME "Quadtree"

#define TWO_ONE 1 // enforce 2:1 refinement ratio

#define I     (point.i - GHOSTS)
#define J     (point.j - GHOSTS)
#define DELTA (1./(1 << point.level))

typedef struct {
  char flags, neighbors;
} Cell;

typedef struct _Quadtree Point;
typedef struct _Quadtree Quadtree;

typedef struct {
  int i, j, level;
} Index;

struct _Quadtree {
  int depth;       /* the maximum depth of the tree */

  Quadtree * back; /* back pointer to the "parent" quadtree */
  char ** m;       /* the grids at each level */
  int i, j, level; /* the current cell index and level */

  Index * index;   /* indices for leaf traversal */
  int nleaves;     /* number of leaves */
  bool dirty;      /* whether indices should be updated */
};

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

#define CELL(m,level,i)  (*((Cell *) &m[level][(i)*(sizeof(Cell) + datasize)]))

/***** Multigrid macros *****/
#define depth()      (((Quadtree *)grid)->depth)
#define aparent(k,l) \
  CELL(point.m, point.level-1, ((point.i+GHOSTS)/2+k)*(_n/2+2*GHOSTS) + \
       (point.j+GHOSTS)/2+l)
#define child(k,l)   \
  CELL(point.m, point.level+1, (2*point.i-GHOSTS+k)*2*(_n + GHOSTS) +	\
       (2*point.j-GHOSTS+l))

/***** Quadtree macros ****/
#define _n (1 << point.level) /* fixme */
#define cell							\
  CELL(point.m, point.level, point.i*(_n + 2*GHOSTS) + point.j)
#define neighbor(k,l)							\
  CELL(point.m, point.level, (point.i + k)*(_n + 2*GHOSTS) + point.j + l)
#define parent             aparent(0,0)
#define alloc_children()						\
  { point.back->dirty = true;						\
    if (point.level == point.depth) alloc_layer(&point); }
#define free_children()    { point.back->dirty = true; }

/***** Data macros *****/
#define data(k,l)							\
  ((double *) &point.m[point.level][((point.i + k)*(_n + 2*GHOSTS) +	\
				     (point.j + l))*(sizeof(Cell) + datasize) \
				    + sizeof(Cell)])
#define field(cell) ((double *)(((char *) &cell) + sizeof(Cell)))
#define _fine(a,k,l) field(child(k,l))[a]
#define _coarse(a,k,l) field(aparent(k,l))[a]

#define POINT_VARIABLES						     \
  VARIABLES							     \
  int level = point.level; NOT_UNUSED(level);			     \
  struct { int x, y; } child = {				     \
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1		     \
  }; NOT_UNUSED(child);

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

#define foreach_boundary_cell(dir)					\
  {									\
    int ig = _ig[dir], jg = _jg[dir];	NOT_UNUSED(ig); NOT_UNUSED(jg);	\
    Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid); \
    int _d = dir; NOT_UNUSED(_d);					\
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;	\
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
	/* do something */
#define end_foreach_boundary_cell()					\
        if (point.level < point.depth) {				\
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

#define foreach_cell() foreach_boundary_cell(nboundary)			\
  POINT_VARIABLES;
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
  }

#define foreach_cell_post(condition)					\
  {									\
    Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid); \
    struct { int l, i, j, stage; } stack[STACKSIZE]; int _s = -1;	\
    _push (0, GHOSTS, GHOSTS, 0); /* the root cell */			\
    while (_s >= 0) {							\
      int stage;							\
      _pop (point.level, point.i, point.j, stage);			\
      switch (stage) {							\
      case 0: {								\
        POINT_VARIABLES;							\
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
        POINT_VARIABLES;							\
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

#define foreach(clause)     { update_cache();				\
  OMP_PARALLEL()							\
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);			\
  Quadtree point = *((Quadtree *)grid); point.back = ((Quadtree *)grid); \
  OMP(omp for schedule(static) clause)					\
  for (int _k = 0; _k < point.nleaves; _k++) {				\
    point.i = point.index[_k].i; point.j = point.index[_k].j;		\
    point.level = point.index[_k].level;				\
    POINT_VARIABLES;
#define end_foreach()         } OMP_END_PARALLEL() }

#define foreach_leaf()            foreach_cell() if (cell.flags & leaf) {
#define end_foreach_leaf()        continue; } end_foreach_cell()

#define foreach_boundary_ghost(dir)        foreach_boundary_cell(dir)	     \
                                             if (cell.flags & leaf) {	     \
					       point.i += ig; point.j += jg; \
					       POINT_VARIABLES;
#define end_foreach_boundary_ghost()   continue; } end_foreach_boundary_cell()

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

Point refine_cell (Point point, scalar * list);

void update_cache (void)
{
  Quadtree * q = grid;
  if (q->dirty) {
    int n = 0;
    foreach_leaf() {
      if (n >= q->nleaves) {
	q->nleaves += 100;
	q->index = realloc (q->index, sizeof (Index)*q->nleaves);
      }
      q->index[n].i = point.i;
      q->index[n].j = point.j;
      q->index[n].level = point.level;
      n++;
    }
    q->nleaves = n;
    q->index = realloc (q->index, sizeof (Index)*q->nleaves);
    q->dirty = false;
  }
}

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
  /* make sure we don't try to access level -1 */
  q->m[0] = NULL; q->m = &(q->m[1]);
  /* initialise the root cell */
  q->m[0] = calloc (_size(0), sizeof (Cell) + datasize);
  CELL(q->m, 0, 2 + 2*GHOSTS).flags |= (leaf | active);
  CELL(q->m, 0, 2 + 2*GHOSTS).neighbors = 1; // only itself as neighbor
  q->index = NULL;
  q->nleaves = 0;
  grid = q;
  while (depth--)
    foreach_leaf()
      point = refine_cell (point, all);
  update_cache();
  init_boundaries (nvar);
  init_events();
}

void free_grid (void)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++)
    free (q->m[l]);
  q->m = &(q->m[-1]);
  free(q->m);
  free(q->index);
  free(q);
  free_boundaries();
}

Point refine_cell (Point point, scalar * list)
{
#if TWO_ONE
  /* refine neighborhood if required */
  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)
      for (int l = 0; l != 2*child.y; l += child.y)
	if (aparent(k,l).flags & leaf) {
	  Point p = point;
	  /* fixme: this should be made
	     independent from the quadtree implementation */
	  p.level = point.level - 1;
	  p.i = (point.i + GHOSTS)/2 + k;
	  p.j = (point.j + GHOSTS)/2 + l;
	  p = refine_cell (p, list);
	  assert (p.m == point.m);
	}
#endif

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
#if 0
      /* bilinear interpolation from coarser level */
      for (scalar v in list)
	fine(v,k,l) = 
	  (9.*v[] + 3.*(v[2*k-1,0] + v[0,2*l-1]) + v[2*k-1,2*l-1])/16.;
#else
      /* linear interpolation from coarser level (conservative) */
      for (scalar v in list)
	fine(v,k,l) = v[] + ((v[1,0] - v[-1,0])*(2*k-1)/8. +
			     (v[0,1] - v[0,-1])*(2*l-1)/8.);
#endif
    }

  return point;
}

void output_cells (FILE * fp);

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
