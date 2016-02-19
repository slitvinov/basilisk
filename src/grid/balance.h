// Dynamic load-balancing

static inline void update_pid (Point point, int pid)
{
  if (cell.pid == pid)
    return;
  if (is_leaf(cell) && cell.neighbors)
    foreach_child()
      cell.pid = pid;
  cell.pid = pid;
  cell.flags &= ~(active|border|ignored);
  if (is_local(cell))
    cell.flags |= active;
}

@def foreach_neighborhood(a, data) {
  int ig = 0; NOT_UNUSED(ig);
#if dimension >= 2
  int jg = 0; NOT_UNUSED(jg);
#endif
#if dimension >= 3
  int kg = 0; NOT_UNUSED(kg);
#endif
  long _len = 0, _alen = (a)->len;
  while (_len < _alen) {
    Point point = *((Point *)((a)->p + _len)); _len += sizeof(Point);
    POINT_VARIABLES;
    Cell * c = (Cell *)((a)->p + _len); NOT_UNUSED(c);
    _len += sizeof(Cell);
    if (c->flags & (data))
      _len += datasize;
@
@def end_foreach_neighborhood()
  }
}
@

static inline void append_neighbor (Point point, Array * a,
				    int pid,
				    const unsigned short next,
				    const unsigned short dnext,
				    const unsigned short lnext,
				    char c, FILE * fp)
{
  const unsigned short flag = 1 << user;
  if (!(cell.flags & (flag|lnext))) {
    array_append (a, &point, sizeof(Point));
    if (cell.flags & dnext)
      array_append (a, &cell, sizeof(Cell) + datasize);
    else
      array_append (a, &cell, sizeof(Cell));
    cell.flags |= flag;
    if (fp)
      fprintf (fp, "%g %g %c %d %d %d %d next\n", x, y, c,
	       (cell.flags & dnext) != 0,
	       (cell.flags & next) != 0,
	       (cell.flags & lnext) != 0,
	       (cell.flags & ignored) != 0);
  }
}

static Array * neighborhood (int pid,
			     const unsigned short next,
			     const unsigned short dnext,
			     const unsigned short lnext,
			     FILE * fp)
{
  Array * a = array_new();
  foreach_cell() {
    if (cell.flags & next) {
      if (level == 0)
	foreach_neighbor()
	  append_neighbor (point, a, pid, next, dnext, lnext, 'a', fp);
      foreach_neighbor(1)
	if (cell.neighbors)
	  foreach_child()
	    append_neighbor (point, a, pid, next, dnext, lnext, 'c', fp);
    }
    else
      continue;
    if (is_leaf(cell))
      continue;
  }
  
  const unsigned short flag = 1 << user;
  foreach_neighborhood(a, dnext) {
    if (is_refined(cell)) {
      c->flags |= next;
      //      int ppid = cell.pid;
      foreach_child()
	if (true/*cell.pid != ppid*/) {
	  // fixme: send only root cells
	  append_neighbor (point, a, pid, next, dnext, lnext, 'e', fp);
	  cell.flags &= ~flag;
	}
    }
    else
      c->flags &= ~next;
    cell.flags &= ~(next|dnext|lnext|flag);
  }

#if 1
  foreach_cell()
    cell.flags &= ~(next|dnext|lnext);
#endif

  return a;
}

static inline void mark_neighbor (Point point,
				  const unsigned short next,
				  const unsigned short dnext)
{
  cell.flags |= next;
  if (is_leaf(cell)) {
    foreach_neighbor()
      cell.flags |= dnext;
    foreach_neighbor(1)
      if (cell.neighbors)
	foreach_child()
	  if (is_rborder(cell) || is_boundary(cell) || is_local(cell))
	    cell.flags |= dnext;
  }
  else { // not leaf
    if (cell.flags & halo)
      cell.flags |= dnext;
    // This is necessary for some schemes i.e. Saint-Venant. The
    // corresponding test case is src/test/bump2Dp.c
    foreach_neighbor(1)
      if (is_leaf(cell))
	cell.flags |= dnext;
  }
}

static bool send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, MOVED_TAG(), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, MOVED_TAG(), MPI_COMM_WORLD, &r[1]);
    return true;
  }
  else
    return false;
}

static bool receive_tree (int from,
			  const unsigned short next, const unsigned short dnext,
			  FILE * fp)
{
  Array a;
  MPI_Recv (&a.len, 1, MPI_LONG, from, MOVED_TAG(),
	    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (a.len > 0) {
    a.p = malloc (a.len);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    MPI_Recv (a.p, a.len, MPI_BYTE, from, MOVED_TAG(),
	      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int n = 0;
    foreach_neighborhood (&a, dnext) {
      if (!(c->flags & ignored)) {
	if (fp)
	  fprintf (fp, "recv %g %g %g %d %d %d %d %d %d %d %d\n",
		   x, y, z, c->pid, cell.pid, from, c->flags & leaf,
		   c->flags & next, c->flags & dnext, c->flags & halo, n);
	update_pid (point, c->pid);
	if (c->flags & dnext)
	  memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
		  datasize);
      }
      if (c->flags & next) {
	if (is_leaf(cell)) {
	  if (fp)
	    fprintf (fp, "refine %g %g %g %d %d %d %d %d %d\n",
		     x, y, z, c->pid, cell.pid, from,
		     c->flags & leaf, cell.flags & leaf, n);
	  refine_cell (point, NULL, 0, NULL, fp);
	}
#if 0
	else
	  foreach_child()
	    if (cell.pid != c->pid) {
	      update_pid (point, c->pid);
	      if (is_local(cell))
		aparent(0).flags |= active;
	    }
#endif
      }
      n++;
    }
    free (a.p);
    return true;
  }
  else
    return false;
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

void check_flags()
{
#if DEBUG_MPI
  foreach_cell()
    foreach_neighbor()
      if (allocated(0))
	for (int i = user; i <= user + 7; i++)
	  assert (!(cell.flags & (1 << i)));
#endif
}

bool balance (double imbalance)
{
  check_flags();

  long nl = 0;
  foreach()
    nl++;
  long nt = nl, nmin = nl, nmax = nl;
  // fixme: do all reductions in one go
  mpi_all_reduce (nt,   MPI_LONG, MPI_SUM);
  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  long ne = max(1, nt/npe());
  //  if (ne < 1000) ne = 1000; // fixme: this does not work

  //  fprintf (stderr, "pid %d nl %ld nt %ld\n", pid(), nl, nt);
  
  if (nmax - nmin <= 1 || nmax - nmin <= ne*imbalance)
    return false;

#if 0  
  {
    char name[80];
    sprintf (name, "before-%d", pid());
    FILE * fp = fopen(name, "w");
    output_cells (fp);

    foreach_cell() {
      fprintf (fp, "pid %g %g %d\n", x, y, cell.pid);
      if (is_leaf(cell))
	continue;
    }

    fclose (fp);
  }
#endif
  
  scalar index[];
  z_indexing (index, true);

  // parallel index restriction
  for (int l = depth(); l >= 0; l--)
    mpi_boundary_restriction_children (mpi_boundary, {index}, l);
    
#if DEBUG_MPI
  char name[80];
  sprintf (name, "bal-%d", pid());
  FILE * fp = fopen (name, "w");
#else
  FILE * fp = NULL;
#endif

  enum {
    prev  = 1 << (user + 1),
    next  = 2*prev,
    dprev = 2*next,
    dnext = 2*dprev,
    lprev = 2*dnext,
    lnext = 2*lprev
  };
  int pid_changed = false;
  foreach_cell_post (!is_leaf(cell)) {

    if (cell.pid == pid() - 1) {
      if (is_leaf(cell) && is_rborder(cell)) {
	foreach_neighbor()
	  if (!is_boundary(cell))
	    cell.flags |= lprev;
      }
      else
	cell.flags |= lprev;
    }
    else if (cell.pid == pid() + 1) {
      if (is_leaf(cell) && is_rborder(cell)) {
	foreach_neighbor()
	  if (!is_boundary(cell))
	    cell.flags |= lnext;
      }
      else
	cell.flags |= lnext;
    }
    
    if (!(cell.flags & ignored)) {
      int pid = balanced_pid (index[], nt);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (cell.pid != pid) {
	if (is_local(cell)) {
	  if (pid == cell.pid - 1)
	    mark_neighbor (point, prev, dprev);
	  else
	    mark_neighbor (point, next, dnext);
	}
	if (fp)
	  fprintf (fp, "changed %g %g %g %d %d %d\n",
		   x, y, z, cell.pid, pid, is_local(cell));
	update_pid (point, pid);
	pid_changed = true;
      }
    }

    if (level > 0) {
      if (cell.flags & prev)
	aparent(0).flags |= prev;
      if (cell.flags & next)
	aparent(0).flags |= next;
    }
  }

  Array * aprev = neighborhood (pid() - 1, prev, dprev, lprev, fp);
  Array * anext = neighborhood (pid() + 1, next, dnext, lnext, fp);

  check_flags();
    
  /* fixme: this should be optimisable using the fact that "send to
     next" must be mutually exclusive with "receive from next" */
  // send mesh to previous/next process
  MPI_Request rprev[2], rnext[2];
  if (pid() > 0         && send_tree (aprev, pid() - 1, rprev))
    quadtree->dirty = true;
  if (pid() < npe() - 1 && send_tree (anext, pid() + 1, rnext))
    quadtree->dirty = true;

  // receive mesh from next/previous process
  if (pid() < npe() - 1 && receive_tree (pid() + 1, prev, dprev, fp))
    quadtree->dirty = true;
  if (pid() > 0 &&         receive_tree (pid() - 1, next, dnext, fp))
    quadtree->dirty = true;  
  
  if (fp)
    fclose (fp);
  
  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);
  
  if (quadtree->dirty) {
#if 1
    // update active cells: fixme: can this be done above
    static const unsigned short refined = 1 << user;
    foreach_cell_post (!is_leaf (cell)) {
      if (!is_leaf(cell) && !is_local(cell)) {
	unsigned short flags = cell.flags & ~active;
	foreach_child()
	  if (is_active(cell))
	    flags |= active, break;
	cell.flags = flags;
      }
      if (cell.flags & refined) // fixme: unused
	cell.flags &= ~refined;
    }
#endif
    flag_border_cells(); // fixme: can this be done above?
    pid_changed = true;
  }

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update();

  return pid_changed;
}
