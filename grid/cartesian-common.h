#include "events.h"

#ifndef is_face_x
# define is_face_x() true
# define is_face_y() true
#endif

#ifndef foreach_boundary_ghost
# define foreach_boundary_ghost(dir)					\
  foreach_boundary(dir,false) {						\
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

#define boundary_flux(...)

void boundary_level (scalar * list, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = _boundary[b][s] (point, s);
}

// Cartesian methods

void (* boundary) (scalar *);

void cartesian_boundary (scalar * list)
{
  boundary_level (list, depth());
}

static double symmetry (Point point, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, scalar s)
{
  return -s[];
}

scalar cartesian_new_scalar (scalar s)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    _boundary[b][s] = symmetry;
  return s;
}

vector cartesian_new_vector (vector v)
{
  /* set default boundary conditions (symmetry) */
  _boundary[top][v.x] = _boundary[bottom][v.x] = symmetry;
  _boundary[right][v.y] = _boundary[left][v.y] = symmetry;
  _boundary[right][v.x] = _boundary[left][v.x] = antisymmetry;
  _boundary[top][v.y] = _boundary[bottom][v.y] = antisymmetry;
  return v;
}

tensor cartesian_new_tensor (tensor t)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    _boundary[b][t.x.x] = _boundary[b][t.y.y] = symmetry;
    _boundary[b][t.x.y] = _boundary[b][t.y.x] = antisymmetry;
  }
  return t;
}

void cartesian_methods()
{
  new_scalar = cartesian_new_scalar;
  new_vector = cartesian_new_vector;
  new_tensor = cartesian_new_tensor;
  boundary   = cartesian_boundary;
}
