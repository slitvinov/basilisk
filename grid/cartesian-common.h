#include "events.h"

#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
  }

#define boundary(...) boundary_level (scalars(__VA_ARGS__), depth())
#define boundary_flux(...)

static bool boundary_var (int b, scalar v, int l)
{
  if (_boundary[b][v]) {
    (* _boundary[b][v]) (l);
    return true;
  }
  return false;
}

void boundary_level (scalar * list, int l)
{
  for (int b = 0; b < nboundary; b++)
    for (scalar s in list) {
      if (!boundary_var (b, s, l))
	foreach_boundary_level (b, l)
	  s[ghost] = s[]; /* default is symmetry */
    }
}

void boundary_uv (scalar u, scalar v)
{
  /* slip walls (symmetry) by default */
  if (!boundary_var (right, u, depth()))
    foreach_boundary (right)
      u[ghost] = 0.;
  if (!boundary_var (right, v, depth()))
    foreach_boundary (right)
      v[ghost] = v[];
  if (!boundary_var (left, u, depth()))
    foreach_boundary (left) {
      u[ghost] = - u[1,0];
      u[] = 0.;
    }
  if (!boundary_var (left, v, depth()))
    foreach_boundary (left)
      v[ghost] = v[];
  if (!boundary_var (top, u, depth()))
    foreach_boundary (top)
      u[ghost] = u[];
  if (!boundary_var (top, v, depth()))
    foreach_boundary (top)
      v[ghost] = 0.;
  if (!boundary_var (bottom, u, depth()))
    foreach_boundary (bottom)
      u[ghost] = u[];
  if (!boundary_var (bottom, v, depth()))
    foreach_boundary (bottom) {
      v[ghost] = - v[0,1];
      v[] = 0.;
    }
}
