#define QUADTREE 1

#include "multigrid-common.h"

// scalar attributes

attribute {
  void (* refine) (Point, scalar);
}

// Quadtree methods

/*
  A cache of refined cells is maintained (if not NULL).
*/
void refine_cell (Point point, scalar * list, int flag, Cache * refined)
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
	  refine_cell (p, list, flag, refined);
	  aparent(k,l).flags |= flag;
	}
    }
#endif

  /* refine */
  int premote = is_remote_leaf(cell) ? remote_leaf : 0;
  cell.flags &= ~(leaf|remote_leaf);

  /* update neighborhood */
  increment_neighbors (point);

  /* for each child: (note that using foreach_child() would be nicer
     but it seems to be significanly slower) */
  int pactive = is_active(cell) ? active : 0;
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      if (dimension == 2)
	assert(!is_active(child(k,l)));
      child(k,l).flags |= (pactive|premote|leaf);
@if _MPI
      if (is_border(cell))
	for (int m = -GHOSTS; m <= GHOSTS; m++)
	  for (int n = -GHOSTS; n <= GHOSTS; n++)
	    if (child(k+m,l+n).pid != cell.pid) {
	      child(k,l).flags |= border;
	      m = n = GHOSTS + 1; // break
	    }
@endif
    }

  if (is_local(cell))
    /* initialise scalars */
    for (scalar s in list)
      s.refine (point, s);

@if _MPI
  if (refined && is_border(cell))
    cache_append (refined, point, 0, 0, cell.flags);
@endif
}

bool coarsen_cell (Point point, scalar * list, CacheLevel * coarsened)
{
#if TWO_ONE
  /* check that neighboring cells are not too fine */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      if (child(k,l).neighbors)
	return false; // cannot coarsen
#endif

  /* restriction */
  for (scalar s in list)
    s.coarsen (point, s);

  /* coarsen */
  cell.flags |= leaf;

  /* update neighborhood */
  decrement_neighbors (point);
  
  /* for each child */
  if (cell.neighbors)
    for (int k = 0; k < 2; k++)
      for (int l = 0; l < 2; l++) {
	child(k,l).flags &= ~(leaf|active);
	child(k,l).pid = cell.pid;
      }

@if _MPI
  if (coarsened && is_border(cell))
    cache_level_append (coarsened, point);
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
    for (int k = -GHOSTS; k < 2 + GHOSTS; k++)
      for (int l = -GHOSTS; l < 2 + GHOSTS; l++)
	if (allocated_child(k,l) &&
	    is_local(child(k,l)) && !is_border(child(k,l)))
	  child(k,l).flags |= border;
  }
@endif

  return true;
}

typedef struct {
  int nc, nf;
} astats;

enum {
  too_coarse = 1 << user,
  too_fine   = 1 << (user + 1)
};

@if _MPI
void mpi_boundary_refine  (void);
void mpi_boundary_coarsen (int);
void mpi_boundary_update  (void);
@else
@ define mpi_boundary_refine()
@ define mpi_boundary_coarsen(a)
@ define mpi_boundary_update()
@endif

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
  scalar * listcm = NULL;

  if (is_constant(cm)) {
    if (p.list == NULL)
      p.list = all;
    restriction (p.slist);
  }
  else {
    if (p.list == NULL) {
      listcm = list_concat (NULL, {cm,fm});
      for (scalar s in all)
	listcm = list_add (listcm, s);
      p.list = listcm;
    }
    scalar * listr = list_concat (p.slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (!is_constant(s) && s.coarsen != no_coarsen)
      listc = list_add (listc, s);

  // refinement
  quadtree->refined.n = 0;
  foreach_cell() {
    if (is_active(cell)) {
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &quadtree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined)
	  // cell has already been refined, skip its children
	  continue;
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	for (int k = 0; k <= 1 && !local; k++)
	  for (int l = 0; l <= 1 && !local; l++)
	    if (is_local(child(k,l)))
	      local = true;
	if (local) {
	  foreach_child()
	    cell.flags &= ~(too_coarse|too_fine);
	  int i = 0;
	  for (scalar s in p.slist) {
	    double max = p.max[i++], sc[4];
	    int c = 0;
	    foreach_child()
	      sc[c++] = s[];
	    s.prolongation (point, s);
	    c = 0;
	    foreach_child() {
	      double e = fabs(sc[c] - s[]);
	      if (e > max) {
		cell.flags &= ~too_fine;
		cell.flags |= too_coarse;
	      }
	      else if (e <= max/1.5 && !(cell.flags & too_coarse))
		cell.flags |= too_fine;
	      s[] = sc[c++];
	    }
	  }
	}
	// cell is too fine, its children cannot be refined
	if (level == p.maxlevel - 1)
	  continue;
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine();
  
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth() - 1; l >= max(p.minlevel, 1); l--) {
    quadtree->coarsened.n = 0;
    foreach_cell() {
      if (is_local (cell)) {
	if (is_leaf (cell))
	  continue;
	else if (level == l) {
	  if (cell.flags & refined)
	    // cell was refined previously, unset the flag
	    cell.flags &= ~(refined|too_fine);
	  else if (cell.flags & too_fine) {
	    bool coarsen = true;
	    foreach_child()
	      if (!(cell.flags & too_fine))
		coarsen = false;
	    if (coarsen && coarsen_cell (point, listc, &quadtree->coarsened))
	      st.nc++;
	    cell.flags &= ~too_fine; // do not coarsen parent
	  }
	  continue;
	}
      } // non-local cell
      else if (level == l || is_leaf(cell))
	continue;
    }
    mpi_boundary_coarsen (l);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf) {
    mpi_boundary_update();
    boundary (p.listb ? p.listb : p.list);
  }
  free (listcm);

  return st;
}

int coarsen_function (int (* func) (Point p), scalar * list)
{
  int nc = 0;
  for (int l = depth() - 1; l >= 0; l--) {
    quadtree->coarsened.n = 0;
    foreach_cell() {
      if (is_local (cell)) { // always true in serial
	if (is_leaf (cell))
	  continue;
	else if (level == l) {
	  if ((*func) (point) &&
	      coarsen_cell (point, list, &quadtree->coarsened))
	    nc++;
	  continue;
	}
      } // non-local cell
      else if (level == l || is_leaf(cell))
	continue;
    }
    mpi_boundary_coarsen (l);
  }
  mpi_boundary_update();
  return nc;
}

#define refine(cond, list) {						\
  quadtree->refined.n = 0;						\
  int refined;								\
  do {									\
    refined = 0;							\
    foreach_leaf()							\
      if (cond)	{							\
        refine_cell (point, list, 0, &quadtree->refined);		\
	refined++;							\
      }									\
  } while (refined);							\
  mpi_boundary_refine();						\
  mpi_boundary_update();						\
}

#define is_refined(cell) ((cell).neighbors && !is_leaf (cell) && \
			  !is_boundary(cell))

static void halo_restriction_flux (vector * list)
{
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
      foreach_halo (prolongation, l) {
	foreach_dimension() {
	  if (is_refined (neighbor(-1,0)))
	    for (vector f in listv)
	      f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
	  if (is_refined (neighbor(1,0)))
	    for (vector f in listv)
	      f.x[1,0] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
        }
      }
    boundary_iterate (halo_restriction_flux, listv);
    free (listv);
  }
}

// Cartesian methods

Point locate (double xp, double yp)
{
  foreach_cell() {
    if (is_active(cell)) {
      Delta /= 2.;
      if (xp < x - Delta || xp > x + Delta || yp < y - Delta || yp > y + Delta)
	continue;
      if (is_leaf(cell) && is_local(cell))
	return point;
    }
    else
      continue;
  }
  Point point = { .level = -1 };
  return point;
}

static scalar quadtree_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  s.refine = s.prolongation = refine_bilinear;
  s.coarsen = coarsen_average;
  return s;
}

foreach_dimension()
static void refine_face_x (Point point, scalar s)
{
  vector v = s.v;
  if (!is_refined(neighbor(-1,0))) {
    double g = (v.x[0,+1] - v.x[0,-1])/8.;
    fine(v.x,0,0) = v.x[] - g;
    fine(v.x,0,1) = v.x[] + g;
  }
  double vx[2];
  if (!is_refined(neighbor(1,0))) {
    double g = (v.x[1,+1] - v.x[1,-1])/8.;
    vx[0] = v.x[1,0] - g;
    vx[1] = v.x[1,0] + g;
    if (neighbor(1,0).neighbors)
      for (int i = 0; i <= 1; i++)
	fine(v.x,2,i) = vx[i];
  }
  else
    for (int i = 0; i <= 1; i++)
      vx[i] = fine(v.x,2,i);
  for (int i = 0; i <= 1; i++)
    fine(v.x,1,i) = (fine(v.x,0,i) + vx[i])/2.;
}

void refine_face (Point point, scalar s)
{
  vector v = s.v;
  foreach_dimension()
    v.x.prolongation (point, v.x);
}

void refine_face_solenoidal (Point point, scalar s)
{
  refine_face (point, s);
  vector v = s.v;
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
  v.y.coarsen = no_coarsen;
  v.x.refine  = refine_face;
  v.y.refine  = no_coarsen;
  v.x.prolongation = refine_face_x;
  v.y.prolongation = refine_face_y;  
  return v;
}

static void quadtree_boundary_level (scalar * list, int l)
{
  if (l < 0)
    l = depth();

  scalar * listdef = NULL, * listc = NULL;
  for (scalar s in list) 
    if (!is_constant (s)) {
      if (s.coarsen == coarsen_average)
	listdef = list_add (listdef, s);
      else if (s.coarsen != no_coarsen)
	listc = list_add (listc, s);
    }

  if (listdef || listc) {
    scalar * list = list_concat (listdef, listc);
    boundary_iterate (halo_restriction, list, l);
    for (int i = l - 1; i >= 0; i--) {
      foreach_halo (restriction, i) {
	for (scalar s in listdef)
	  s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
	for (scalar s in listc)
	  s.coarsen (point, s);
      }
      boundary_iterate (halo_restriction, list, i);
    }
    free (list);
    free (listdef);
    free (listc);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  for (scalar s in list)
    if (!is_constant (s) && s.refine != no_coarsen) {
      if (s.face)
	listf = vectors_add (listf, s.v);
      else
	listr = list_add (listr, s);
    }

  if (listr || listf) {
    boundary_iterate (halo_prolongation, list, 0, l);
    for (int i = 0; i < l; i++) {
      foreach_halo (prolongation, i) {
	for (scalar s in listr)
          s.prolongation (point, s);
	for (vector v in listf)
	  foreach_dimension()
	    v.x.prolongation (point, v.x);
      }
      boundary_iterate (halo_prolongation, list, i + 1, l);
    }
    free (listr);
    free (listf);
  }
}

void quadtree_methods()
{
  multigrid_methods();
  init_scalar      = quadtree_init_scalar;
  init_face_vector = quadtree_init_face_vector;
  boundary_level   = quadtree_boundary_level;
  boundary_flux    = halo_restriction_flux;
}
