#include <math.h>

#define foreach_fine_to_coarse()           foreach_cell_post(!(cell.flags & leaf))
#define end_foreach_fine_to_coarse()       end_foreach_cell_post()

#define foreach_level(l)                   foreach_cell() { \
                                             if (level == l || cell.flags & leaf) {
#define end_foreach_level()                  continue; } } end_foreach_cell()

#define foreach_boundary(dir)              foreach_boundary_cell(dir)	\
                                             if (cell.flags & leaf) {	\
                                               QUADTREE_VARIABLES;	\
					       VARIABLES;
#define end_foreach_boundary()               continue; } end_foreach_boundary_cell()

#define foreach_boundary_level(dir,l)      foreach_boundary_cell(dir)               \
                                             QUADTREE_VARIABLES;	            \
                                             if (level == l || cell.flags & leaf) { \
					       VARIABLES;
#define end_foreach_boundary_level()         continue; } end_foreach_boundary_cell()

#undef boundary_ghost
#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
    int _in = -_ig[d], _jn = -_jg[d];					\
    foreach_halo() if (neighbor(_in,_jn).flags & leaf) { x; }		\
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
    if (cell.flags & leaf)
      return point;
  }
  Point point = {-1, NULL, NULL, -1, -1, -1}; // not found
  return point;
}

#include "multigrid-common.h"

bool coarsen_cell (Point point)
{
  QUADTREE_VARIABLES;

#if TWO_ONE
  /* check that neighboring cells are not too fine */
  for (int k = -1; k < 3; k++)
    for (int l = -1; l < 3; l++)
      if ((child(k,l).flags & active) && !(child(k,l).flags & leaf))
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
      if (cell.flags & leaf)
	continue;
      else if (level == l) {
	if ((*func) (point) && coarsen_cell (point))
	  nc++;
	continue;
      }
    }
  return nc;
}

int coarsen_wavelet (scalar w, double max)
{
  int nc = 0;
  for (int l = depth() - 1; l >= 0; l--)
    foreach_cell() {
      if (cell.flags & leaf)
	continue;
      else if (level == l) {
	double error = 0.;
	for (int k = 0; k < 2; k++)
	  for (int l = 0; l < 2; l++) {
	    double e = fabs(fine(w,k,l));
	    if (e > error)
	      error = e;
	  }
	if (error < max && coarsen_cell(point))
	  nc++;
	/* propagate the error to coarser levels */
	w[] = fabs(w[]) + error;	
	continue;
      }
    }
  return nc;
}

int refine_function (scalar start, scalar end,
		     int (* func) (Point p, void * data), 
		     void * data)
{
  int nf = 0;
  foreach_leaf()
    if ((*func) (point, data)) {
      point = refine_cell (point, start, end);
      nf++;
    }
  return nf;
}

int refine_wavelet (scalar start, scalar end,
		    scalar w, double max)
{
  int nf = 0;
  foreach_leaf()
    /* fixme: w[] should be explicitly defined */
    if (w[] != undefined && fabs(w[]) >= max) {
      point = refine_cell (point, start, end);
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++)
	  fine(w,k,l) = undefined;
      nf++;
    }
  return nf;
}

int flag_halo_cells ()
{
  int nh = 0;

  /* reset old halos first */
  foreach_cell() {
    if (!(cell.flags & halo))
      continue;
    else {
      cell.flags &= ~halo;
#if TRASH
      if (!(cell.flags & active))
	for (scalar v = 0; v < nvar; v++)
	  val(v,0,0) = undefined;
#endif
    }
  }

  /* from the bottom up */
  foreach_cell_post (cell.neighbors > 0 || (cell.flags & active)) {
    if (!(cell.flags & active)) {
      /* inactive and neighbors > 0 => this is a halo cell */
      cell.flags |= halo;
      /* propagate to parent */
      parent.flags |= halo;
      nh++;
    }
    else if ((cell.flags & halo) && level > 0)
      /* propagate to parent */
      parent.flags |= halo;
  }

  return nh;
}

#define foreach_halo() foreach_cell() { \
  if (!(cell.flags & halo))		      \
    continue;				      \
  else if (!(cell.flags & active)) {
#define end_foreach_halo()  }} end_foreach_cell();

/* breadth-first traversal of halos from coarse to fine */
#define foreach_halo_coarse_fine(depth1)    {				\
  int _depth = depth1 < 0 ? depth() : depth1;				\
  for (int _l = 0; _l <= _depth; _l++)                                  \
    foreach_cell() {							\
      if (!(cell.flags & halo))                                         \
      	continue; /* no more halos, skip the rest of this branch */     \
      if (level == _l) {                                                \
	if (!(cell.flags & active))
#define end_foreach_halo_coarse_fine()                                  \
	continue; /* already at level l, skip the deeper branches */    \
      }                                                                 \
    } end_foreach_cell();				                \
}

void update_halo (int depth, scalar start, scalar end)
{
  foreach_halo_coarse_fine (depth)
    /* bilinear interpolation from coarser level */
    for (scalar v = start; v <= end; v++)
      val(v,0,0) = (9.*coarse(v,0,0) + 
		    3.*(coarse(v,childx,0) + coarse(v,0,childy)) + 
		    coarse(v,childx,childy))/16.;
}

void update_halo_u_v (int depth, scalar u, scalar v)
{
  foreach_halo_coarse_fine (depth) {
    /* linear interpolation from coarser level */
    if (childx < 0)
      /* conservative interpolation */
      u[] = coarse(u,0,0) + (coarse(u,0,1) - coarse(u,0,-1))*childy/8.;
    else
      u[] = (3.*coarse(u,0,0) + coarse(u,0,childy) + 
	     3.*coarse(u,1,0) + coarse(u,1,childy))/8.;
    if (childy < 0)
      /* conservative interpolation */
      v[] = coarse(v,0,0) + (coarse(v,1,0) - coarse(v,-1,0))*childx/8.;
    else
      v[] = (3.*coarse(v,0,0) + coarse(v,childx,0) + 
	     3.*coarse(v,0,1) + coarse(v,childx,1))/8.;
  }
}

#undef boundary
void boundary (scalar p)
{
  boundary_level (p, depth());
  restriction (p, p);
  update_halo (-1, p, p);
}

#undef boundary_flux
void boundary_flux (vector f)
{
  restriction_flux (f);
  foreach_cell() {
    if (!(cell.flags & halo))
      continue;
    else if (cell.flags & leaf) {
      if (child(0,0).flags & halo) {
	if (child(0,1).flags & halo)
	  f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
	if (child(1,0).flags & halo)
	  f.y[] = (fine(f.y,0,0) + fine(f.y,1,0))/2.;
      }
      continue;
    }
  }
}