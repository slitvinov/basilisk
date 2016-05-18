#define QUADTREE 1 // OBSOLETE: for backward-compatibility
#define TREE 1

#include "multigrid-common.h"

// scalar attributes

attribute {
  void (* refine) (Point, scalar);
}

// Tree methods

/*
  A cache of refined cells is maintained (if not NULL).
*/
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{
  int nr = 0;
#if TWO_ONE
  /* refine neighborhood if required */
  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)
    #if dimension > 1
      for (int l = 0; l != 2*child.y; l += child.y)
    #endif
      #if dimension > 2
	for (int m = 0; m != 2*child.z; m += child.z)
      #endif
	  if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
	    Point p = point;
	    /* fixme: this should be made
	       independent from the tree implementation */
	    p.level = point.level - 1;
	    p.i = (point.i + GHOSTS)/2 + k;
            #if dimension > 1
	      p.j = (point.j + GHOSTS)/2 + l;
            #endif
            #if dimension > 2
	      p.k = (point.k + GHOSTS)/2 + m;
            #endif
	      nr += refine_cell (p, list, flag, refined);
	    aparent(k,l,m).flags |= flag;
	  }
#endif

  /* refine */
  cell.flags &= ~leaf;

  /* update neighborhood */
  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
  foreach_child()
    cell.flags |= cflag;
    
  /* initialise scalars */
  for (scalar s in list)
    if (is_local(cell) || s.face)
      s.refine (point, s);

@if _MPI
  if (is_border(cell)) {
    foreach_child() {
      bool bord = false;
      foreach_neighbor() {
	if (!is_local(cell) || (level > 0 && !is_local(aparent(0))))
	  bord = true, break;
	if (is_refined(cell))
	  foreach_child()
	    if (!is_local(cell))
	      bord = true, break;
	if (bord)
	  break;
      }
      if (bord)
	cell.flags |= border;
    }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
@endif
  return nr;
}

bool coarsen_cell (Point point, scalar * list)
{
#if TWO_ONE
  /* check that neighboring cells are not too fine.
     check that children are not different boundaries */
  int pid = cell.pid;
  foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; // cannot coarsen
#endif

  /* restriction */
  for (scalar s in list)
    s.coarsen (point, s);

  /* coarsen */
  cell.flags |= leaf;

  /* update neighborhood */
  decrement_neighbors (point);
  
@if _MPI
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
    foreach_neighbor(1)
      if (cell.neighbors)
	foreach_child()
	  if (allocated(0) && is_local(cell) && !is_border(cell))
	    cell.flags |= border;
  }
@endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{
#if TWO_ONE
  /* recursively coarsen children cells */
  foreach_child()
    if (cell.neighbors)
      foreach_neighbor(1)
	if (is_refined (cell))
	  coarsen_cell_recursive (point, list);
#endif
  assert (coarsen_cell (point, list));
}

void mpi_boundary_refine  (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update  (scalar *);

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
  if (p.minlevel < 1)
    p.minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &tree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined) {
	  // cell has already been refined, skip its children
	  cell.flags &= ~too_coarse;
	  continue;
	}
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	if (!local)
	  foreach_child()
	    if (is_local(cell))
	      local = true, break;
	if (local) {
	  int i = 0;
	  static const int just_fine = 1 << (user + 3);
	  for (scalar s in p.slist) {
	    double max = p.max[i++], sc[1 << dimension];
	    int c = 0;
	    foreach_child()
	      sc[c++] = s[];
	    s.prolongation (point, s);
	    c = 0;
	    foreach_child() {
	      double e = fabs(sc[c] - s[]);
	      if (e > max && level < p.maxlevel) {
		cell.flags &= ~too_fine;
		cell.flags |= too_coarse;
	      }
	      else if ((e <= max/1.5 || level > p.maxlevel) &&
		       !(cell.flags & (too_coarse|just_fine))) {
		if (level >= p.minlevel)
		  cell.flags |= too_fine;
	      }
	      else if (!(cell.flags & too_coarse)) {
		cell.flags &= ~too_fine;
		cell.flags |= just_fine;
	      }
	      s[] = sc[c++];
	    }
	  }
	  foreach_child() {
	    cell.flags &= ~just_fine;
	    if (!is_leaf(cell)) {
	      cell.flags &= ~too_coarse;
	      if (level >= p.maxlevel)
		cell.flags |= too_fine;
	    }
	    else if (!is_active(cell))
	      cell.flags &= ~too_coarse;
	  }
	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);
  
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= p.minlevel; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
	if (level == l) {
	  if (!is_leaf(cell)) {
	    if (cell.flags & refined)
	      // cell was refined previously, unset the flag
	      cell.flags &= ~(refined|too_fine);
	    else if (cell.flags & too_fine) {
	      if (is_local(cell) && coarsen_cell (point, listc))
		st.nc++;
	      cell.flags &= ~too_fine; // do not coarsen parent
	    }
	  }
	  if (cell.flags & too_fine)
	    cell.flags &= ~too_fine;
	  else if (aparent(0).flags & too_fine)
	    aparent(0).flags &= ~too_fine;
	  continue;
	}
	else if (is_leaf(cell))
	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);
  free (listcm);
  
  return st;
}

#define refine(cond, list) do {			                        \
  int refined;								\
  do {									\
    refined = 0;							\
    tree->refined.n = 0;						\
    foreach_leaf()							\
      if (cond) {							\
	refine_cell (point, list, 0, &tree->refined);			\
	refined++;							\
	continue;							\
      }									\
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);				\
    if (refined) {							\
      mpi_boundary_refine (list);					\
      mpi_boundary_update (list);					\
    }									\
  } while (refined);							\
} while(0)

static void refine_level (int depth)
{
  refine (level < depth, NULL);
}

#define unrefine(cond, list) do {					\
  static const int too_fine = 1 << user;			        \
  foreach_cell() {							\
    if (is_leaf(cell))							\
      continue;								\
    if (is_local(cell) && (cond))					\
      cell.flags |= too_fine;						\
  }									\
  for (int _l = depth(); _l >= 0; _l--) {				\
    foreach_cell() {							\
      if (is_leaf(cell))						\
	continue;							\
      if (level == _l) {						\
	if (is_local(cell) && (cell.flags & too_fine)) {		\
	  coarsen_cell (point, list);					\
	  cell.flags &= ~too_fine;					\
	}								\
	continue;							\
      }									\
    }									\
    mpi_boundary_coarsen (_l, too_fine);				\
  }									\
  mpi_boundary_update (list);						\
} while (0)

static void halo_flux (vector * list)
{
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
      foreach_halo (prolongation, l)
	foreach_dimension() {
#if dimension == 1
	  if (is_refined (neighbor(-1)))
	    for (vector f in listv)
	      f.x[] = fine(f.x,0);
	  if (is_refined (neighbor(1)))
	    for (vector f in listv)
	      f.x[1] = fine(f.x,2);
#elif dimension == 2
	  if (is_refined (neighbor(-1)))
	    for (vector f in listv)
	      f.x[] = (fine(f.x,0,0) + fine(f.x,0,1))/2.;
	  if (is_refined (neighbor(1)))
	    for (vector f in listv)
	      f.x[1] = (fine(f.x,2,0) + fine(f.x,2,1))/2.;
#else // dimension == 3
	  if (is_refined (neighbor(-1)))
	    for (vector f in listv)
	      f.x[] = (fine(f.x,0,0,0) + fine(f.x,0,1,0) +
		       fine(f.x,0,0,1) + fine(f.x,0,1,1))/4.;
	  if (is_refined (neighbor(1)))
	    for (vector f in listv)
	      f.x[1] = (fine(f.x,2,0,0) + fine(f.x,2,1,0) +
			fine(f.x,2,0,1) + fine(f.x,2,1,1))/4.;
#endif
      }
    free (listv);
  }
}

// Cartesian methods

static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  s.refine = s.prolongation;
  return s;
}

foreach_dimension()
static void refine_face_x (Point point, scalar s)
{
  assert (dimension > 1);
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {
    double g1 = (v.x[0,+1] - v.x[0,-1])/8.;
    double g2 = (v.x[0,0,+1] - v.x[0,0,-1])/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(v.x,0,j,k) = v.x[] + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {
    double g1 = (v.x[1,+1] - v.x[1,-1])/8.;
    double g2 = (v.x[1,0,+1] - v.x[1,0,-1])/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(v.x,2,j,k) = v.x[1] + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (v.x[0,+1] + v.x[1,+1] - v.x[0,-1] - v.x[1,-1])/16.;
    double g2 = (v.x[0,0,+1] + v.x[1,0,+1] - v.x[0,0,-1] - v.x[1,0,-1])/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(v.x,1,j,k) = (v.x[] + v.x[1])/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
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
#if dimension > 1
  if (is_local(cell)) {
    // local projection, see section 3.3 of Popinet, JCP, 2009
    vector v = s.v;
    double d[1 << dimension], p[1 << dimension];
    int i = 0;
    foreach_child() {
      d[i] = 0.;
      foreach_dimension()
	d[i] += v.x[1] - v.x[];
      i++;
    }
#if dimension == 2
    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1) += p[1] - p[0];
    fine(v.x,1,0) += p[3] - p[2];
    fine(v.y,0,1) += p[0] - p[2];
    fine(v.y,1,1) += p[1] - p[3];
#elif dimension == 3
    static double m[7][7] = {
      {7./12,5./24,3./8,5./24,3./8,1./4,1./3},
      {5./24,7./12,3./8,5./24,1./4,3./8,1./3},
      {3./8,3./8,3./4,1./4,3./8,3./8,1./2},
      {5./24,5./24,1./4,7./12,3./8,3./8,1./3},
      {3./8,1./4,3./8,3./8,3./4,3./8,1./2},
      {1./4,3./8,3./8,3./8,3./8,3./4,1./2},
      {1./3,1./3,1./2,1./3,1./2,1./2,5./6}
    };
    p[0] = 0.;
    for (int i = 0; i < 7; i++) {
      p[i + 1] = 0.;
      for (int j = 0; j < 7; j++)
	p[i + 1] += m[i][j]*d[j+1];
    }
    for (int k = 0; k <= 1; k++) {
      fine(v.x,1,0,k) += p[4+k] - p[0+k];
      fine(v.x,1,1,k) += p[6+k] - p[2+k];
      fine(v.y,0,1,k) += p[2+k] - p[0+k];
      fine(v.y,1,1,k) += p[6+k] - p[4+k];
    }
    fine(v.z,0,0,1) += p[1] - p[0];
    fine(v.z,0,1,1) += p[3] - p[2];
    fine(v.z,1,0,1) += p[5] - p[4];
    fine(v.z,1,1,1) += p[7] - p[6];
#endif // dimension == 3
  }
#endif // dimension > 1
}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  foreach_dimension()
    v.x.coarsen = v.x.refine = no_coarsen;
  v.x.coarsen = coarsen_face;
  v.x.refine  = refine_face;
  foreach_dimension()
    v.x.prolongation = refine_face_x;
  return v;
}

static void tree_boundary_level (scalar * list, int l)
{
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    boundary_iterate (level, list, depth);
    return;
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  for (scalar s in list) 
    if (!is_constant (s)) {
      if (s.coarsen == coarsen_average) {
	listdef = list_add (listdef, s);
	list2 = list_add (list2, s);
      }
      else if (s.coarsen != no_coarsen) {
	listc = list_add (listc, s);
	if (s.face)
	  foreach_dimension()
	    list2 = list_add (list2, s.v.x);
	else
	  list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    boundary_iterate (restriction, list2, depth);
    for (int l = depth - 1; l >= 0; l--) {
      foreach_coarse_level(l) {
	for (scalar s in listdef)
	  coarsen_average (point, s);
	for (scalar s in listc)
	  s.coarsen (point, s);
      }
      boundary_iterate (restriction, list2, l);
    }
    free (listdef);
    free (listc);
    free (list2);
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
    boundary_iterate (level, list, 0);
    for (int i = 0; i < depth; i++) {
      foreach_halo (prolongation, i) {
	for (scalar s in listr)
          s.prolongation (point, s);
	for (vector v in listf)
	  foreach_dimension()
	    v.x.prolongation (point, v.x);
      }
      boundary_iterate (level, list, i + 1);
    }
    free (listr);
    free (listf);
  }
}

double treex (Point point) {
  if (level == 0)
    return 0;
#if dimension == 2
  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;
#else
  assert (false); // not implemented
  double i = 0;
#endif
  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) {
  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
  foreach_cell()
    if (cell.neighbors)
      foreach_child()
	if (is_local(cell))
	  fprintf (fp, "%g %g\n%g %g\n\n",
		   treex(parent), treey(parent), treex(point), treey(point));
}

void tree_check()
{
  // checks the consistency of the tree

  long nleaves = 0, nactive = 0;
  foreach_cell_all() {
    if (is_leaf(cell)) {
      assert (cell.pid >= 0); // boundaries cannot be leaves
      nleaves++;
    }
    if (is_local(cell))
      assert (is_active(cell) || is_prolongation(cell));
    if (is_active(cell))
      nactive++;
    // check number of refined neighbors
    int neighbors = 0;
    foreach_neighbor(1)
      if (allocated(0) && is_refined(cell))
	neighbors++;
    assert (cell.neighbors == neighbors);
    
    // checks that prolongation cells do not have children
    if (!cell.neighbors)
      assert (!allocated_child(0));
  }

  // checks that all active cells are reachable
  long reachable = 0;
  foreach_cell() {
    if (is_active(cell))
      reachable++;
    else
      continue;
  }
  assert (nactive == reachable);

  // checks that all leaf cells are reachable
  reachable = 0;
  foreach_cell()
    if (is_leaf(cell)) {
      reachable++;
      continue;
    }
  assert (nleaves == reachable);
}

static void tree_restriction (scalar * list) {
  if (tree_is_full())
    multigrid_restriction (list);
  // otherwise the restriction has already been done by tree_boundary_level()
}

void tree_methods()
{
  multigrid_methods();
  init_scalar      = tree_init_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level   = tree_boundary_level;
  boundary_flux    = halo_flux;
  restriction      = tree_restriction;
}