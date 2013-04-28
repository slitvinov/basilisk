#include "events.h"

#ifndef foreach_boundary_ghost
# define foreach_boundary_ghost(dir)					\
  foreach_boundary(dir) {						\
    point.i += ig; point.j += jg;					\
    ig = -ig; jg = -jg;							\
    POINT_VARIABLES;
# define end_foreach_boundary_ghost()					\
    ig = -ig; jg = -jg;							\
  } end_foreach_boundary()
#endif

#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
  }

#define boundary(...) boundary_level (scalars(__VA_ARGS__), depth())
#define boundary_flux(...)

typedef double (* BoundaryFunc) (Point point, scalar s);

void boundary_level (scalar * list, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = (* ((BoundaryFunc)_boundary[b][s])) (point, s);
}

static double symmetry (Point point, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, scalar s)
{
  return -s[];
}

scalar new_scalar (scalar s)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    _boundary[b][s] = symmetry;
  return s;
}

vector new_vector (vector v)
{
  /* set default boundary conditions (symmetry) */
  _boundary[top][v.x] = _boundary[bottom][v.x] = symmetry;
  _boundary[right][v.y] = _boundary[left][v.y] = symmetry;
  _boundary[right][v.x] = _boundary[left][v.x] = antisymmetry;
  _boundary[top][v.y] = _boundary[bottom][v.y] = antisymmetry;
  return v;
}

tensor new_tensor (tensor t)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    _boundary[b][t.x.x] = _boundary[b][t.y.y] = symmetry;
    _boundary[b][t.x.y] = _boundary[b][t.y.x] = antisymmetry;
  }
  return t;
}
