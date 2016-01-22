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
#define _BACK   (2*point.k - GHOSTS)
#define _FRONT  (_BACK + 1)

#if 0
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
#endif

#define STACKSIZE 20
#if dimension == 1
#define _push(b,c,d,e,f)			  \
  { _s++; stack[_s].l = b;			  \
    stack[_s].i = c;				  \
    stack[_s].stage = f; }
#define _pop()					  \
  { point.level = stack[_s].l;			  \
    point.i = stack[_s].i;			  \
    stage = stack[_s].stage; _s--; }
#elif dimension == 2
#define _push(b,c,d,e,f)			  \
  { _s++; stack[_s].l = b;			  \
    stack[_s].i = c; stack[_s].j = d;		  \
    stack[_s].stage = f; }
#define _pop()					  \
  { point.level = stack[_s].l;			  \
    point.i = stack[_s].i; point.j = stack[_s].j; \
    stage = stack[_s].stage; _s--; }
#else // dimension == 3
#define _push(b,c,d,e,f)						\
  { _s++; stack[_s].l = b;						\
    stack[_s].i = c; stack[_s].j = d; stack[_s].k = e;			\
    stack[_s].stage = f; }
#define _pop()								\
  { point.level = stack[_s].l;						\
    point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; \
    stage = stack[_s].stage; _s--; }
#endif

@def foreach_cell()
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 1
    Point point = {GHOSTS,0};
    struct { int l, i, stage; } stack[STACKSIZE];
#elif dimension == 2
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    Point point = {GHOSTS,GHOSTS,GHOSTS,0};
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, GHOSTS, GHOSTS, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0))
	continue;
      switch (stage) {
      case 0: {
	POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell()
        if (point.level < quadtree->depth) {
	  _push (point.level, point.i, point.j, point.k, 1);
          _push (point.level + 1, _LEFT, _TOP, _FRONT, 0);
        }
        break;
      }
#if dimension == 1
      case 1: _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0); break;
#else // dimension >= 2
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	      _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	      _push (point.level + 1, _LEFT,  _BOTTOM, _FRONT, 0); break;
#endif
#if dimension == 2
      case 3: _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0); break;
#else // dimension >= 3
      case 3: _push (point.level, point.i, point.j, point.k, 4);
	      _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0); break;
      case 4: _push (point.level, point.i, point.j, point.k, 5);
	      _push (point.level + 1, _LEFT, _TOP, _BACK, 0); break;
      case 5: _push (point.level, point.i, point.j, point.k, 6);
	      _push (point.level + 1, _RIGHT, _TOP, _BACK, 0); break;
      case 6: _push (point.level, point.i, point.j, point.k, 7);
	      _push (point.level + 1, _LEFT, _BOTTOM, _BACK, 0); break;
      case 7: _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0); break;
#endif
      }
    }
  }
@

@def foreach_cell_post(condition)
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
#if dimension == 1
    Point point = {GHOSTS,0};
    struct { int l, i, stage; } stack[STACKSIZE];
#elif dimension == 2
    Point point = {GHOSTS,GHOSTS,0};
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    Point point = {GHOSTS,GHOSTS,GHOSTS,0};
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, GHOSTS, GHOSTS, GHOSTS, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
	if (point.level == quadtree->depth) {
	  _push (point.level, point.i, point.j, point.k, 8);
	}
	else {
	  _push (point.level, point.i, point.j, point.k, 1);
	  if (condition)
	    _push (point.level + 1, _LEFT, _TOP, _FRONT, 0);
	}
	break;
      }
      case 1:
	_push (point.level, point.i, point.j, point.k, 2);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0);
	break;
#if dimension >= 2
      case 2:
	_push (point.level, point.i, point.j, point.k, 3);
	if (condition)
	  _push (point.level + 1, _LEFT,  _BOTTOM, _FRONT, 0);
	break;
      case 3:
	_push (point.level, point.i, point.j, point.k, 4);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0);
	break;
#endif
#if dimension >= 3
      case 4:
	_push (point.level, point.i, point.j, point.k, 5);
	if (condition)
	  _push (point.level + 1, _LEFT, _TOP, _BACK, 0);
	break;
      case 5:
	_push (point.level, point.i, point.j, point.k, 6);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _BACK, 0);
	break;
      case 6:
	_push (point.level, point.i, point.j, point.k, 7);
	if (condition)
	  _push (point.level + 1, _LEFT,  _BOTTOM, _BACK, 0);
	break;
      case 7:
	_push (point.level, point.i, point.j, point.k, 8);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0);
	break;	
#endif
      default: {
        POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell_post()
      }
      }
    }
  }
@
