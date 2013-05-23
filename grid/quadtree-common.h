#define QUADTREE 1

#ifndef foreach
# define foreach foreach_leaf
# define end_foreach end_foreach_leaf
#endif

#ifndef foreach_fine_to_coarse
# define foreach_fine_to_coarse()      foreach_cell_post(!is_leaf(cell))
# define end_foreach_fine_to_coarse()  end_foreach_cell_post()
#endif

#ifndef foreach_level
# define foreach_level(l)						\
  foreach_cell() {							\
    if (!is_active (cell))						\
      continue;								\
    else if (level == l) {
# define end_foreach_level()    continue; } } end_foreach_cell()
#endif

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf(l)   foreach_cell() { \
                                      if (level == l || is_leaf(cell)) {
# define end_foreach_level_or_leaf()  continue; } } end_foreach_cell()
#endif

#define foreach_halo()     foreach_halo_coarse_to_fine(-1)
#define end_foreach_halo() end_foreach_halo_coarse_to_fine()

#define foreach_boundary(dir,corners)				\
  foreach_boundary_cell(dir,corners) if (is_leaf (cell)) {
#define end_foreach_boundary()  continue; } end_foreach_boundary_cell()

#define foreach_boundary_level(dir,l,corners)	\
  foreach_boundary_cell(dir,corners)		\
    if (level == l || is_leaf(cell)) {
#define end_foreach_boundary_level()				\
      corners(); /* we need this otherwise we'd skip corners */	\
      continue;							\
    }								\
  end_foreach_boundary_cell()

#define is_face_x() !is_refined(neighbor(-1,0))
#define is_face_y() !is_refined(neighbor(0,-1))

#include "multigrid-common.h"

// Quadtree methods

int coarsen_function (int (* func) (Point p), scalar * list)
{
  int nc = 0;
  for (int l = depth() - 1; l >= 0; l--)
    foreach_cell() {
      if (is_leaf (cell))
	continue;
      else if (level == l) {
	if ((*func) (point) && coarsen_cell (point, list))
	  nc++;
	continue;
      }
    }
  return nc;
}

int coarsen_wavelet (scalar w, double max, int minlevel, scalar * list)
{
  int nc = 0;
  scalar * listc = NULL;
  for (scalar s in list)
    if (s.coarsen != refine_none)
      listc = list_append (listc, s);
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
	// propagate the error to coarser levels
	double wc = fabs(w[]) + error;
	if (error < max && level >= minlevel && coarsen_cell (point, listc))
	  nc++;
	w[] = wc;
	continue;
      }
    }
  free (listc);
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

static void huge (Point point, scalar w)
{
  // prevents coarsening of newly created cells (when combined
  // with coarsen_wavelet())
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      fine(w,k,l) = HUGE;  
}

int refine_wavelet (scalar w, double max, int maxlevel, scalar * list)
{
  // overload the refine method for w
  void * f = w.refine;
  w.refine = huge;
  scalar * list1 = NULL;
  for (scalar s in list)
    if (s.refine != refine_none)
      list1 = list_append (list1, s);
  if (!list_lookup (list1, w))
    // add w to the list of variables to refine
    list1 = list_append (list1, w);
  // refine
  int nf = 0;
  foreach_leaf()
    if (level < maxlevel && w[] != HUGE && fabs(w[]) >= max) {
      point = refine_cell (point, list1);
      nf++;
    }
  free (list1);
  // restore refine method
  w.refine = f;
  return nf;
}

void halo_restriction (scalar * def, scalar * list)
{
  foreach_halo_fine_to_coarse () {
    for (scalar s in def)
      s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
    for (scalar s in list)
      s.coarsen (point, s);
  }
}

void halo_restriction_flux (vector * list)
{
  foreach_halo_fine_to_coarse()
    foreach_dimension() {
      for (vector f in list)
        f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
      if (is_leaf (neighbor(1,0)))
	for (vector f in list)
	  f.x[1,0] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
    }
}

void halo_prolongation (int depth, scalar * list)
{
  foreach_halo_coarse_to_fine (depth) {
    for (scalar s in list)
      if (s.gradient) { // linear interpolation (e.g. with limiting)
	struct { double x, y; } g;
	foreach_dimension()
	  g.x = s.gradient (coarse(s,-1,0), 
			    coarse(s,0,0), 
			    coarse(s,1,0));
	s[] = coarse(s,0,0) + (g.x*child.x + g.y*child.y)/4.;
      }
      else // fixme: should this use s.refine?
	/* bilinear interpolation from coarser level */
	s[] = (9.*coarse(s,0,0) + 
	       3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
	       coarse(s,child.x,child.y))/16.;	
  }
}

void halo_prolongation_u_v (int depth, scalar u, scalar v)
{
  foreach_halo_coarse_to_fine (depth) {
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

// Multigrid methods

void quadtree_boundary_restriction (scalar * list)
{
  // traverse the boundaries of all coarse levels (i.e. not the leaves)
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_cell (b, true) {
      if (is_leaf (cell))
	continue;
      else
	for (scalar s in list)
	  s[ghost] = s.boundary[b] (point, s);
    }
}

// Cartesian methods

#undef boundary_ghost
#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
    int _in = -_ig[d], _jn = -_jg[d];					\
    foreach_halo() if (is_leaf(_neighbor(_in,_jn))) {			\
      ig = _in; jg = _jn; VARIABLES; x; }				\
    end_foreach_halo();							\
  }

#undef boundary_flux
#define boundary_flux halo_restriction_flux

void quadtree_boundary (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL;
  for (scalar s in list)
    if (s.coarsen == coarsen_average)
      listdef = list_append (listdef, s);
    else if (s.coarsen != refine_none)
      listc = list_append (listc, s);

  if (listdef || listc) {
    halo_restriction (listdef, listc);
    for (int b = 0; b < nboundary; b++)
      foreach_boundary_cell (b, true) {
	if (is_active (cell)) {
	  if (cell.neighbors > 0) {
	    for (scalar s in listdef)
	      s[ghost] = s.boundary[b] (point, s);
	    for (scalar s in listc)
	      s[ghost] = s.boundary[b] (point, s);
	  }
	}
	else
	  continue;
      }
    free (listdef);
    free (listc);
  }

  scalar * listr = NULL;
  for (scalar s in list)
    if (s.refine != refine_none)
      listr = list_append (listr, s);

  if (listr) {
    halo_prolongation (-1, listr);
    for (int b = 0; b < nboundary; b++)
      foreach_boundary_cell (b, true)
	if (!is_active (cell)) {
	  if (cell.neighbors > 0) {
	    for (scalar s in listr)
	      s[ghost] = s.boundary[b] (point, s);
	  }
	  else
	    continue;
	}
    free (listr);
  }
}

Point locate (double xp, double yp)
{
  foreach_cell () {
    delta /= 2.;
    if (xp < x - delta || xp > x + delta || yp < y - delta || yp > y + delta)
      continue;
    if (is_leaf (cell))
      return point;
  }
  Point point = { level: -1 };
  return point;
}

scalar quadtree_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  s.refine = refine_linear;
  s.coarsen = coarsen_average;
  return s;
}

void quadtree_methods()
{
  multigrid_methods();
  init_scalar          = quadtree_init_scalar;
  boundary             = quadtree_boundary;
  boundary_restriction = quadtree_boundary_restriction;
}
