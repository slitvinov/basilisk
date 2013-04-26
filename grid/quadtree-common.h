#define foreach_fine_to_coarse()		\
  foreach_cell_post(!is_leaf(cell))
#define end_foreach_fine_to_coarse()		\
  end_foreach_cell_post()
#define foreach_level(l)		        \
  foreach_cell() {				\
  if (level == l || is_leaf(cell)) {
#define end_foreach_level()			\
  continue; } } end_foreach_cell()

#define foreach_boundary(dir)		        \
  foreach_boundary_cell(dir)			\
  if (is_leaf(cell)) {				\
    POINT_VARIABLES;
#define end_foreach_boundary()			\
  continue; } end_foreach_boundary_cell()

#define foreach_boundary_level(dir,l)		\
  foreach_boundary_cell(dir)			\
    POINT_VARIABLES;				\
    if (level == l || is_leaf(cell)) {		\
      POINT_VARIABLES;
#define end_foreach_boundary_level()            \
    continue; } end_foreach_boundary_cell()

#include "multigrid-common.h"

#undef boundary_ghost
#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
    int _in = -_ig[d], _jn = -_jg[d];					\
    foreach_halo() if (is_leaf(neighbor(_in,_jn))) { x; }		\
    end_foreach_halo();							\
  }

Point locate (double xp, double yp)
{
  foreach_cell () {
    double delta = DELTA;
    double x = (point.i - GHOSTS + 0.5)*delta - 0.5;
    double y = (point.j - GHOSTS + 0.5)*delta - 0.5;
    delta /= 2.;
    if (xp < x - delta || xp > x + delta || yp < y - delta || yp > y + delta)
      continue;
    if (is_leaf (cell))
      return point;
  }
  Point point = {-1, NULL, NULL, -1, -1, -1}; // not found
  return point;
}

bool coarsen_cell (Point point)
{
#if TWO_ONE
  /* check that neighboring cells are not too fine */
  for (int k = -1; k < 3; k++)
    for (int l = -1; l < 3; l++)
      if (is_active (child(k,l)) && !is_leaf (child(k,l)))
	return false; // cannot coarsen
#endif

  /* coarsen */
  free_children();
  cell.flags |= leaf;
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      child(k,l).flags &= ~(leaf|active);
#if TRASH
      /* trash the data just to make sure it's never touched */
      for (scalar v = 0; v < nvar; v++)
	fine(v,k,l) = undefined;
#endif
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors--;
    }
  return true;
}

int coarsen_function (int (* func) (Point p))
{
  int nc = 0;
  for (int l = depth() - 1; l >= 0; l--)
    foreach_cell() {
      if (is_leaf (cell))
	continue;
      else if (level == l) {
	if ((*func) (point) && coarsen_cell (point))
	  nc++;
	continue;
      }
    }
  return nc;
}

int coarsen_wavelet (scalar w, double max, int minlevel)
{
  int nc = 0;
  for (int l = depth() - 1; l >= 0; l--)
    foreach_cell() {
      if (is_leaf (cell))
	continue;
      else if (level == l) {
	double error = 0.;
	for (int k = 0; k < 2; k++)
	  for (int l = 0; l < 2; l++) {
	    double e = fabs(fine(w,k,l));
	    if (e > error)
	      error = e;
	  }
	if (error < max && level >= minlevel && coarsen_cell(point))
	  nc++;
	/* propagate the error to coarser levels */
	w[] = fabs(w[]) + error;	
	continue;
      }
    }
  return nc;
}

int refine_function (int (* func) (Point p, void * data), 
		     void * data,
		     scalar * list)
{
  int nf = 0;
  foreach_leaf()
    if ((*func) (point, data)) {
      point = refine_cell (point, list);
      nf++;
    }
  return nf;
}

int refine_wavelet (scalar w, double max, int maxlevel, scalar * list)
{
  int nf = 0;
  foreach_leaf()
    /* fixme: w[] should be explicitly defined */
    if (w[] != undefined && fabs(w[]) >= max && level < maxlevel) {
      point = refine_cell (point, list);
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++)
	  fine(w,k,l) = undefined;
      nf++;
    }
  return nf;
}

void update_halo (int depth, scalar * list)
{
  foreach_halo_coarse_fine (depth)
    /* bilinear interpolation from coarser level */
    for (scalar s in list)
      s[] = (9.*coarse(s,0,0) + 
	     3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
	     coarse(s,child.x,child.y))/16.;
}

void update_halo_u_v (int depth, scalar u, scalar v)
{
  foreach_halo_coarse_fine (depth) {
    /* linear interpolation from coarser level */
    if (child.x < 0)
      /* conservative interpolation */
      u[] = coarse(u,0,0) + (coarse(u,0,1) - coarse(u,0,-1))*child.y/8.;
    else
      u[] = (3.*coarse(u,0,0) + coarse(u,0,child.y) + 
	     3.*coarse(u,1,0) + coarse(u,1,child.y))/8.;
    if (child.y < 0)
      /* conservative interpolation */
      v[] = coarse(v,0,0) + (coarse(v,1,0) - coarse(v,-1,0))*child.x/8.;
    else
      v[] = (3.*coarse(v,0,0) + coarse(v,child.x,0) + 
	     3.*coarse(v,0,1) + coarse(v,child.x,1))/8.;
  }
}

#undef boundary
#define boundary(...) boundary_a(scalars(__VA_ARGS__))
void boundary_a (scalar * list)
{
  boundary_level (list, depth());
  restriction (list);
  update_halo (-1, list);
}

#undef boundary_flux
#define boundary_flux(...) boundary_flux_a(vectors(__VA_ARGS__))
void boundary_flux_a (vector * list)
{
  restriction_flux (list);
  foreach() {
    if (is_active(neighbor(-1,0)) && !is_leaf (neighbor(-1,0))) {
      for (vector v in list)
	v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
    }
    if (is_active(neighbor(0,-1)) && !is_leaf (neighbor(0,-1))) {
      for (vector v in list)
	v.y[] = (fine(v.y,0,0) + fine(v.y,1,0))/2.;
    }
  }
}
