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
@ def foreach_level_or_leaf(l)
  foreach_cell() {
    if (level == l || is_leaf(cell)) {
@
@ define end_foreach_level_or_leaf()  continue; } } end_foreach_cell()
@endif

@define foreach_halo()     foreach_halo_coarse_to_fine(depth())
@define end_foreach_halo() end_foreach_halo_coarse_to_fine()

@def foreach_boundary_halo(dir)
  foreach_boundary_cell (dir, false)
    if (!is_active(cell)) {
      if (cell.neighbors > 0) {
	point.i += ig; point.j += jg;
	ig = -ig; jg = -jg;
	POINT_VARIABLES;
@
@def end_foreach_boundary_halo()
        ig = -ig; jg = -jg;
        point.i -= ig; point.j -= jg;
      }
      else
	continue;
    }
  end_foreach_boundary_cell()
@

@def foreach_boundary(dir,corners)
  foreach_boundary_cell(dir,corners) if (is_leaf (cell)) { @
@def end_foreach_boundary()  continue; } end_foreach_boundary_cell() @

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
  foreach_halo() if (allocated(_in,_jn) && is_leaf(neighbor(_in,_jn))) {
    ig = _in; jg = _jn; VARIABLES;
@
@define end_foreach_boundary_ghost_halo() } end_foreach_halo(); }

@def foreach_boundary_fine_to_coarse(dir)
  foreach_boundary_cell_post (dir, !is_leaf (cell)) {
    point.i += ig; point.j += jg;
    ig = -ig; jg = -jg;
    POINT_VARIABLES;
@
@def end_foreach_boundary_fine_to_coarse()
    ig = -ig; jg = -jg;
  } end_foreach_boundary_cell_post()
@

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
	  aparent(k,l).flags |= refined;
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
    // if (is_leaf (neighbor(-1,0)))
      for (vector f in list)
	f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
      if (is_leaf (neighbor(1,0)))
	for (vector f in list)
	  f.x[1,0] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
    }
}

void halo_prolongation (int depth, scalar * list)
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

void quadtree_boundary_centered (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL;
  for (scalar s in list) {
    if (s.coarsen == coarsen_average)
      listdef = list_append (listdef, s);
    else if (s.coarsen != refine_none)
      listc = list_append (listc, s);
  }

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
    halo_prolongation (depth(), listr);
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

void quadtree_boundary_tangent (vector * list)
{
  /* we need to define the normal ghost values on coarse (right/top) boundaries
     (this is not done by boundary_normal()) */
  foreach_boundary_fine_to_coarse(right)
    for (vector u in list)
      u.x[] = (fine(u.x,0,0) + fine(u.x,0,1))/2.;
  foreach_boundary_fine_to_coarse(top)
    for (vector u in list)
      u.y[] = (fine(u.y,0,0) + fine(u.y,1,0))/2.;

  cartesian_boundary_tangent (list);

  // interior halos
  foreach_halo()
    foreach_dimension()
      // fixme: this does not work for static quadtrees 
      // because allocated() is always true
      if (allocated(-1,0) && !is_leaf(neighbor(-1,0))) {
	if (child.x < 0) {
	  for (vector u in list)
	    u.x[] = (3.*coarse(u.x,0,0) + coarse(u.x,0,child.y))/4.;
	}
	else
	  for (vector u in list)
	    u.x[] = (3.*coarse(u.x,0,0) + coarse(u.x,0,child.y) +
		     3.*coarse(u.x,1,0) + coarse(u.x,1,child.y))/8.;
      }

  // we need to do the same for the (right/top) boundary halos
  foreach_boundary_halo (right)
    for (vector u in list)
      u.x[] = (3.*coarse(u.x,0,0) + coarse(u.x,0,child.y))/4.;
  foreach_boundary_halo (top)
    for (vector u in list)
      u.y[] = (3.*coarse(u.y,0,0) + coarse(u.y,child.x,0))/4.;
}

Point locate (double xp, double yp)
{
  foreach_cell () {
    Delta /= 2.;
    if (xp < x - Delta || xp > x + Delta || yp < y - Delta || yp > y + Delta)
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
  assert (d[0] + d[1] + d[2] + d[3] < 1e-3);
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
  v = multigrid_init_face_vector (v, name);
  v.x.refine  = refine_face;
  v.y.refine  = coarsen_none;
  return v;
}

void quadtree_methods()
{
  multigrid_methods();
  init_scalar           = quadtree_init_scalar;
  boundary_centered     = quadtree_boundary_centered;
  boundary_normal       = halo_restriction_flux;
  boundary_tangent      = quadtree_boundary_tangent;
  boundary_restriction  = quadtree_boundary_restriction;
  init_face_vector = quadtree_init_face_vector;
}
