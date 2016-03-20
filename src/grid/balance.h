// Dynamic load-balancing

@if TRASH
@ define is_indexed(v) (cell.pid < 0 || (!isnan(v) && v >= 0))
@else
@ define is_indexed(v) (cell.pid < 0 || v >= 0)
@endif

Array * neighborhood (scalar index, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
  foreach_cell() {
    // root cells
    bool root = false;
    if ((!is_local(cell) || index[] != nextpid) && is_refined(cell)) {
      foreach_child()
	if (is_local(cell) && index[] == nextpid)
	  root = true, break;
      if (root && cell.pid != nextpid) {
	if (fp)
	  fprintf (fp, "%g %g %g %d root\n", x, y, index[], cell.pid);
	foreach_neighbor()
	  if (cell.pid != nextpid && is_indexed(index[]))
	    cell.flags |= sent;
      }
    }
    // children
    if ((is_local(cell) && index[] == nextpid) || root) {
      if (fp)
	fprintf (fp, "%g %g %g %d nextpid\n", x, y, index[], cell.pid);
      foreach_neighbor(1)
	if (cell.neighbors && cell.pid != nextpid)
	  foreach_child()
	    if (cell.pid != nextpid && is_indexed(index[]))
	      cell.flags |= sent;
    }
    if (is_leaf(cell))
      continue;
  }

  return tree (sizeof(Cell) + datasize);
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

static bool receive_tree (int from, scalar index, FILE * fp)
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
    const unsigned short next = 1 << (user + 1);
    foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
	      datasize);
      assert (cell.pid < 0 || index[] >= 0);
      if (fp)
	fprintf (fp, "%g %g %g %d %d %d %d recv\n",
		 x, y, index[], cell.pid,
		 c->flags & leaf,
		 cell.flags & leaf, from);
      if ((c->flags & next) && (c->flags & leaf) && !is_leaf(cell)) {
	if (fp)
	  fprintf (fp, "%g %g %g %d blarg\n", x, y, index[], cell.pid);
	if (cell.neighbors)
	  assert (coarsen_cell (point, NULL));
	cell.flags |= leaf;
      }
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

static void check_flags()
{
#if DEBUG_MPI
  foreach_cell()
    foreach_neighbor()
      if (allocated(0))
	for (int i = user; i <= user + 7; i++)
	  assert (!(cell.flags & (1 << i)));
#endif
}

static void check_s()
{
#if DEBUG_MPI
  extern double t;
  if (t <= 0.1)
    return;
  scalar s = {1};
  foreach()
    foreach_neighbor()
      if (!is_boundary(cell))
	assert (!isnan(s[]));
#endif
}

trace
bool balance (double imbalance)
{
  check_flags();
  check_s();
  
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
  
  scalar index[];
  z_indexing (index, true);
  
#if DEBUG_MPI
  char name[80];
  sprintf (name, "bal-%d", pid());
  FILE * fp = fopen (name, "w");
#else
  FILE * fp = NULL;
#endif

  // compute new pid, stored in index[]
  foreach_cell() {
    if (is_local(cell)) {
      int pid = balanced_pid (index[], nt);
      index[] = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (fp)
	fprintf (fp, "%g %g %g %d index\n", x, y, index[], cell.pid);
    }
    else
      index[] = -1;
  }
  for (int l = 0; l <= depth(); l++)
    boundary_iterate (restriction, {index}, l);

  Array * anext = neighborhood (index, pid() + 1, fp);
  Array * aprev = neighborhood (index, pid() - 1, fp);

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
  if (pid() < npe() - 1 && receive_tree (pid() + 1, index, fp))
    quadtree->dirty = true;
  if (pid() > 0 &&         receive_tree (pid() - 1, index, fp))
    quadtree->dirty = true;

  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  foreach_cell_all() {
    if (fp)
      fprintf (fp, "%g %g %g %d %d nindex\n", x, y, index[], cell.pid,
	       is_leaf(cell));
    if (level < depth() && !cell.neighbors &&
	point.i > 0 && point.i <= (1 << level) + 2 &&
	point.j > 0 && point.j <= (1 << level) + 2 &&
	allocated_child(0)) {
      if (fp)
	fprintf (fp, "%g %g %g %d freechildren\n", x, y, index[], cell.pid);
      free_children (point);
    }
  }

  // set new pids
  int pid_changed = false;
  foreach_cell() {
    if (is_indexed(index[]) && cell.pid != index[]) {
      if (fp)
	fprintf (fp, "%g %g %g %d %d %d new\n", x, y, index[], cell.pid,
		 is_leaf(cell), cell.neighbors);
      cell.pid = index[];
      cell.flags &= ~(active|border);
      if (is_local(cell))
	cell.flags |= active;
      if (is_leaf(cell) && cell.neighbors) {
	int pid = cell.pid;
	foreach_child()
	  cell.pid = pid;
      }
      pid_changed = true;
    }
    if (is_leaf(cell))
      continue;
  }

  if (quadtree->dirty || pid_changed) {
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
  
  if (fp)
    fclose (fp);

  check_s();
  
  return pid_changed;
}
