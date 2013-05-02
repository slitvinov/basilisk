#include "cartesian-common.h"

// Multigrid methods

void (* boundary_level)       (scalar *, int);
void (* boundary_restriction) (scalar *);

void restriction (scalar * list)
{
  foreach_fine_to_coarse()
    for (scalar s in list)
      s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
}

void wavelet (scalar v, scalar w)
{
  restriction (v);
  boundary_restriction (v);
  foreach_fine_to_coarse() {
    /* difference between fine value and bilinearly-interpolated
       coarse value */
    fine(w,0,0) = fine(v,0,0) - 
      (9.*v[] + 3.*(v[-1,0] + v[0,-1]) + v[-1,-1])/16.;
    fine(w,0,1) = fine(v,0,1) - 
      (9.*v[] + 3.*(v[-1,0] + v[0,+1]) + v[-1,+1])/16.;
    fine(w,1,0) = fine(v,1,0) - 
      (9.*v[] + 3.*(v[+1,0] + v[0,-1]) + v[+1,-1])/16.;
    fine(w,1,1) = fine(v,1,1) - 
      (9.*v[] + 3.*(v[+1,0] + v[0,+1]) + v[+1,+1])/16.;
  }
  /* root cell */
  foreach_level(0) w[] = 0.;
}

void refine_bilinear (Point point, scalar v)
{
  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      /* bilinear interpolation from coarser level */
      fine(v,k,l) = (9.*v[] + 
		     3.*(v[2*k-1,0] + v[0,2*l-1]) + 
		     v[2*k-1,2*l-1])/16.;
}

void refine_linear (Point point, scalar v)
{
  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      /* linear interpolation from coarser level (conservative) */
      fine(v,k,l) = v[] + ((v[1,0] - v[-1,0])*(2*k-1)/8. +
			   (v[0,1] - v[0,-1])*(2*l-1)/8.);
}

void refine_reset (Point point, scalar v)
{
  /* foreach_child() */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      fine(v,k,l) = 0.;
}

void multigrid_boundary_level (scalar * list, int l)
{
  // traverse the boundaries of a given level
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = _boundary[b][s] (point, s);
}

void multigrid_boundary_restriction (scalar * list)
{
  // traverse the boundaries of all coarse levels (i.e. not the leaves)
  for (int l = 0; l < depth(); l++)
    boundary_level (list, l);
}

void multigrid_methods()
{
  cartesian_methods();
  boundary_level = multigrid_boundary_level;
  boundary_restriction = multigrid_boundary_restriction;
}
