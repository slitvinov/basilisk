#define QUADTREE 1

#include "multigrid-common.h"

// scalar attributes

attribute {
  void   (* refine)       (Point, scalar);
  void   (* coarsen)      (Point, scalar);
};

// Quadtree methods

Point refine_cell (Point point, scalar * list)
{
#if TWO_ONE
  /* refine neighborhood if required */
  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x) {
#if dimension == 1
      int l = 0;
#else
      for (int l = 0; l != 2*child.y; l += child.y)
#endif
	if (aparent(k,l).flags & leaf) {
	  Point p = point;
	  /* fixme: this should be made
	     independent from the quadtree implementation */
	  p.level = point.level - 1;
	  p.i = (point.i + GHOSTS)/2 + k;
	  p.j = (point.j + GHOSTS)/2 + l;
	  p = refine_cell (p, list);
	  assert (p.m == point.m);
	  aparent(k,l).flags |= refined;
	}
    }
#endif

  /* refine */
  cell.flags &= ~leaf;

  /* update neighborhood */
  increment_neighbors (&point);

  /* for each child: (note that using foreach_child() would be nicer
     but it seems to be significanly slower) */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      if (dimension == 2)
	assert(!(child(k,l).flags & active));
      child(k,l).flags |= (active|leaf);
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

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      child(k,l).flags &= ~(leaf|active);

  /* update neighborhood */
  decrement_neighbors (&point);

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
  scalar * listb; // fields which need BCs (default list)
};

astats adapt_wavelet (struct Adapt p)
{
  if (p.list == NULL)
    p.list = all;

  restriction (p.slist);

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (s.coarsen != none)
      listc = list_add (listc, s);

  // refinement
  foreach_cell() {
    if (is_leaf (cell)) {
      if (level < p.maxlevel) {
	int i = 0, refine = false;
	for (scalar s in p.slist) {
	  double w = s.prolongation ?
	    /* difference between fine value and its prolongation */
	    s[] - s.prolongation (point, s) :
	    /* difference between fine value and bilinearly-interpolated
	       coarse value */
	    s[] - (9.*coarse(s,0,0) + 
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
      continue;
    }
    // !is_leaf (cell)
    else if (cell.flags & refined)
      // cell has already been refined, skip its children
      continue;
  }

  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth() - 1; l >= max(p.minlevel, 1); l--)
    foreach_cell() {
      if (is_leaf (cell))
	continue;
      else if (level == l) {
	if (cell.flags & refined) {
	  // cell was refined previously, unset the flag and skip its children
	  cell.flags &= ~refined;
	  continue;
	}
	int i = 0, coarsen = true;
	for (scalar s in p.slist) {
	  double error;
	  if (s.prolongation) {
	    /* difference between fine value and its prolongation */
	    error = fabs (s[] - s.prolongation (point, s));
	    foreach_child() {
	      double e = fabs(s[] - s.prolongation (point, s));
	      if (e > error)
		error = e;
	    }
	  }
	  else {
	    /* difference between fine value and bilinearly-interpolated
	       coarse value */
	    error = fabs (s[] - (9.*coarse(s,0,0) + 
				 3.*(coarse(s,child.x,0) + 
				     coarse(s,0,child.y)) + 
				 coarse(s,child.x,child.y))/16.);
	    for (int k = 0; k < 2; k++)
	      for (int l = 0; l < 2; l++) {
		double e = fabs(fine(s,k,l) - (9.*s[] + 
					       3.*(s[2*k-1,0] + s[0,2*l-1]) + 
					       s[2*k-1,2*l-1])/16.);
		if (e > error)
		  error = e;
	      }
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
    boundary (p.listb ? p.listb : p.list);

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

int refine_function (int (* func) (Point point, void * data), 
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

static void halo_restriction (scalar * def, scalar * listc, int l)
{
  scalar * list = list_concat (def, listc);
  boundary_iterate (halo_restriction, list, l);
  for (l--; l >= 0; l--) {
    foreach_halo (restriction, l) {
      for (scalar s in def)
	s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
      for (scalar s in listc)
	s.coarsen (point, s);
    }
    boundary_iterate (halo_restriction, list, l);
  }
  free (list);
}

static void halo_restriction_flux (vector * list)
{
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv)
    for (int l = depth() - 1; l >= 0; l--)
      foreach_halo (restriction, l)
	foreach_dimension() {
          // if (is_leaf (neighbor(-1,0)))
          for (vector f in listv)
	    f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
	  if (is_leaf (neighbor(1,0)))
	    for (vector f in listv)
	      f.x[1,0] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
        }

  free (listv);
}

static void halo_prolongation (scalar * list, int depth)
{
  for (scalar s in list)
    if (s.gradient)
      s.refine = refine_linear; // fixme: this should be done automatically

  boundary_iterate (halo_prolongation, list, 0, 0);
  for (int l = 0; l < depth; l++) {
    foreach_halo (prolongation, l)
      for (scalar s in list)
	s.refine (point, s);
    boundary_iterate (halo_prolongation, list, l + 1, depth);
  }
}

static void halo_prolongation (scalar * list, int depth)
{
  foreach_halo_coarse_to_fine (depth)
    for (scalar s in list) {
      if (s.gradient) { // linear interpolation (e.g. with limiting)
	double sc = coarse(s,0,0);
	s[] = sc;
	foreach_dimension()
	  s[] += s.gradient (coarse(s,-1,0), sc, coarse(s,1,0))*child.x/4.;
      }
      else if (s.prolongation) // variable-specific prolongation (e.g. VOF)
	s[] = s.prolongation (point, s);
      else
	/* default is bilinear interpolation from coarser level */
	s[] = (9.*coarse(s,0,0) + 
	       3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
	       coarse(s,child.x,child.y))/16.;	
    }
  boundary_halo_prolongation (list, depth);
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

Point locate (double xp, double yp)
{
  foreach_cell() 
    if (is_active(cell)) {
      Delta /= 2.;
      if (xp < x - Delta || xp > x + Delta || yp < y - Delta || yp > y + Delta)
	continue;
      if (is_leaf (cell))
	return point;
    }
  Point point = { .level = -1 };
  return point;
}

static scalar quadtree_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  s.refine  = refine_bilinear;
  s.coarsen = coarsen_average;
  return s;
}

static void refine_face (Point point, scalar s)
{
  vector v = s.v;
  foreach_dimension() {
    if (is_leaf(neighbor(-1,0)) || is_ghost(neighbor(-1,0))) {
      double g = (v.x[0,+1] - v.x[0,-1])/8.;
      fine(v.x,0,0) = v.x[] - g;
      fine(v.x,0,1) = v.x[] + g;
    }
    if (is_leaf(neighbor(+1,0)) || is_ghost(neighbor(+1,0))) {
      double g = (v.x[1,+1] - v.x[1,-1])/8.;
      fine(v.x,2,0) = v.x[1,0] - g;
      fine(v.x,2,1) = v.x[1,0] + g;
    }
    fine(v.x,1,0) = (fine(v.x,0,0) + fine(v.x,2,0))/2.;
    fine(v.x,1,1) = (fine(v.x,0,1) + fine(v.x,2,1))/2.;
  }

  // local projection, see section 3.3 of Popinet, JCP, 2009
  double d[4], p[4];
  d[0] = fine(v.x,1,1) - fine(v.x,0,1) + fine(v.y,0,2) - fine(v.y,0,1);
  d[1] = fine(v.x,2,1) - fine(v.x,1,1) + fine(v.y,1,2) - fine(v.y,1,1);
  d[2] = fine(v.x,1,0) - fine(v.x,0,0) + fine(v.y,0,1) - fine(v.y,0,0);
  d[3] = fine(v.x,2,0) - fine(v.x,1,0) + fine(v.y,1,1) - fine(v.y,1,0);

  p[0] = 0.;
  p[1] = (3.*d[1] + d[2])/4. + d[3]/2.;
  p[2] = (d[1] + 3.*d[2])/4. + d[3]/2.;
  p[3] = (d[1] + d[2])/2. + d[3];
  fine(v.x,1,1) += p[1] - p[0];
  fine(v.x,1,0) += p[3] - p[2];
  fine(v.y,0,1) += p[0] - p[2];
  fine(v.y,1,1) += p[1] - p[3];
}

vector quadtree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  v.x.coarsen = coarsen_face;
  v.y.coarsen = none;
  v.x.refine = refine_face;
  v.y.refine = none;
  return v;
}

static void quadtree_boundary_level (scalar * list, int l)
{
  if (l < 0)
    l = depth();

  scalar * listdef = NULL, * listc = NULL;
  for (scalar s in list) {
    if (s.coarsen == coarsen_average)
      listdef = list_add (listdef, s);
    else if (s.coarsen != none)
      listc = list_add (listc, s);
  }

  if (listdef || listc) {
    halo_restriction (listdef, listc, l);
    free (listdef);
    free (listc);
  }

  scalar * listr = NULL;
  for (scalar s in list)
    if (s.refine != none)
      listr = list_add (listr, s);

  if (listr) {
    halo_prolongation (listr, l);
    free (listr);
  }
}

void quadtree_methods()
{
  multigrid_methods();
  init_scalar          = quadtree_init_scalar;
  init_face_vector     = quadtree_init_face_vector;
  boundary_level       = quadtree_boundary_level;
  boundary_flux        = halo_restriction_flux;
}
