#define QUADTREE 1

@ifndef foreach
@ define foreach foreach_leaf
@ define end_foreach end_foreach_leaf
@endif

@ifndef foreach_fine_to_coarse
@ define foreach_fine_to_coarse()      foreach_cell_post(!is_leaf(cell))
@ define end_foreach_fine_to_coarse()  end_foreach_cell_post()
@endif

@ifndef foreach_level
@ def foreach_level(l)
  foreach_cell() {
    if (!is_active (cell))
      continue;
    else if (level == l) {
@
@ define end_foreach_level()    continue; } } end_foreach_cell()
@endif

@ifndef foreach_level_or_leaf
@ define foreach_level_or_leaf(l)   foreach_cell() { \
                                      if (level == l || is_leaf(cell)) {
@ define end_foreach_level_or_leaf()  continue; } } end_foreach_cell()
@endif

@define foreach_halo()     foreach_halo_coarse_to_fine(-1)
@define end_foreach_halo() end_foreach_halo_coarse_to_fine()

@define foreach_boundary(dir,corners)				\
  foreach_boundary_cell(dir,corners) if (is_leaf (cell)) {
@define end_foreach_boundary()  continue; } end_foreach_boundary_cell()

@def foreach_boundary_level(dir,l,corners)
  foreach_boundary_cell(dir,corners)
    if (level == l || is_leaf(cell)) {
@
@def end_foreach_boundary_level()
      corners(); /* we need this otherwise we'd skip corners */
      continue;
    }
  end_foreach_boundary_cell()
@
  
@def foreach_boundary_ghost_halo(d) {
  int _in = -_ig[d], _jn = -_jg[d];
  foreach_halo() if (_ALLOCATED(_in,_jn) && is_leaf(_neighbor(_in,_jn))) {
    ig = _in; jg = _jn; VARIABLES;
@
@define end_foreach_boundary_ghost_halo() \
  } end_foreach_halo(); }

@define is_face_x() !is_refined(neighbor(-1,0))
@define is_face_y() !is_refined(neighbor(0,-1))

#include "multigrid-common.h"

// Quadtree methods

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

  /* refine */
  cell.flags &= ~leaf;
  /* update neighborhood */
  for (int o = -GHOSTS; o <= GHOSTS; o++)
    for (int p = -GHOSTS; p <= GHOSTS; p++)
      neighbor(o,p).neighbors--;

  /* for each child: (note that using foreach_child() would be nicer
     but it seems to be significanly slower) */
  alloc_children (&point);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      assert(!(child(k,l).flags & active));
      child(k,l).flags |= (active | leaf);
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors++;
    }

  /* initialise scalars */
  for (scalar s in list)
    s.refine (point, s);

  return point;
}

bool coarsen_cell (Point point, scalar * list)
{
#if TWO_ONE
  /* check that neighboring cells are not too fine */
  for (int k = -1; k < 3; k++)
    for (int l = -1; l < 3; l++)
      if (is_active (child(k,l)) && !is_leaf (child(k,l)))
	return false; // cannot coarsen
#endif

  /* restriction */
  for (scalar s in list)
    s.coarsen (point, s);

  /* coarsen */
  cell.flags |= leaf;
  /* update neighborhood */
  for (int o = -GHOSTS; o <= GHOSTS; o++)
    for (int p = -GHOSTS; p <= GHOSTS; p++)
      neighbor(o,p).neighbors++;

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      child(k,l).flags &= ~(leaf|active);
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors--;
    }

  free_children (&point);
  return true;
}

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist; // list of scalars
  double * max;   // tolerance for each scalar
  int maxlevel;   // maximum level of refinement
  int minlevel;   // minimum level of refinement (default 1)
  scalar * list;  // list of fields to update (default all)
};

astats adapt_wavelet (struct Adapt p)
{
  if (p.list == NULL)
    p.list = all;

  restriction (p.slist);
  boundary_restriction (p.slist);

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (s.coarsen != refine_none)
      listc = list_append (listc, s);

  // refinement
  foreach_leaf()
    if (level < p.maxlevel) {
      int i = 0, refine = false;
      for (scalar s in p.slist) {
	/* difference between fine value and bilinearly-interpolated
	   coarse value */
	double w = s[] - (9.*coarse(s,0,0) + 
			  3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
			  coarse(s,child.x,child.y))/16.;
	if (fabs(w) > p.max[i++]) {
	  refine = true; // error is too large
	  break; // no need to check other fields
	}
      }
      if (refine) {
	point = refine_cell (point, p.list);
	st.nf++;
      }
    }

  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth() - 1; l >= max(p.minlevel, 1); l--)
    foreach_cell() {
      if (is_leaf (cell))
	continue;
      else if (level == l) {
	int i = 0, coarsen = true;
	/* difference between fine value and bilinearly-interpolated
	   coarse value */
	for (scalar s in p.slist) {
	  double error =
	    fabs (s[] - (9.*coarse(s,0,0) + 
			 3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
			 coarse(s,child.x,child.y))/16.);
	  for (int k = 0; k < 2; k++)
	    for (int l = 0; l < 2; l++) {
	      double e = fabs(fine(s,k,l) - (9.*s[] + 
					     3.*(s[2*k-1,0] + s[0,2*l-1]) + 
					     s[2*k-1,2*l-1])/16.);
	      if (e > error)
		error = e;
	    }
	  if (error > p.max[i++]/1.5) {
	    coarsen = false; // error is too large
	    break; // no need to check other fields
	  }
	}
	if (coarsen && coarsen_cell (point, listc))
	  st.nc++;
	continue;
      }
    }
  free (listc);

  if (st.nc || st.nf)
    boundary (p.list);

  return st;
}

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
	double sc = coarse(s,0,0);
	s[] = sc;
	foreach_dimension()
	  s[] += s.gradient (coarse(s,-1,0), sc, coarse(s,1,0))*child.x/4.;
      }
      else
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
  Point point = { .level = -1 };
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
