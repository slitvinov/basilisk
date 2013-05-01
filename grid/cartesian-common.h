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

#define boundary(...) boundary_level (scalars(__VA_ARGS__), depth())
#define boundary_flux(...)

void boundary_level (scalar * list, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = (*boundary[b][s]) (point, s);
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
    boundary[b][s] = symmetry;
  return s;
}

vector cartesian_new_vector (vector v)
{
  /* set default boundary conditions (symmetry) */
  boundary[top][v.x] = boundary[bottom][v.x] = symmetry;
  boundary[right][v.y] = boundary[left][v.y] = symmetry;
  boundary[right][v.x] = boundary[left][v.x] = antisymmetry;
  boundary[top][v.y] = boundary[bottom][v.y] = antisymmetry;
  return v;
}

tensor cartesian_new_tensor (tensor t)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    boundary[b][t.x.x] = boundary[b][t.y.y] = symmetry;
    boundary[b][t.x.y] = boundary[b][t.y.x] = antisymmetry;
  }
  return t;
}

void cartesian_methods()
{
  new_scalar = cartesian_new_scalar;
  new_vector = cartesian_new_vector;
  new_tensor = cartesian_new_tensor;
}
