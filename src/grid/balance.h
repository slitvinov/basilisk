// Dynamic load-balancing

typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;

@define NEWPID() ((NewPid *)&val(newpid,0,0,0))

/* fixme: !isnan(v) is only useful because of problems in
   box_boundaries() on level 0. */
@define is_newpided() (!isnan(val(newpid,0,0,0)) && NEWPID()->pid > 0)

Array * tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);  
  Array * a = array_new();

  foreach_cell_post_all (true)
    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0).flags |= next;

  bool empty = true;
  foreach_cell_all() {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && NEWPID()->leaf)
	assert (is_leaf(cell));
      if (is_refined(cell)) {
	/* check that we are not too refined i.e. that none of our
	   children are prolongations */
	bool prolo = false;
	foreach_child()
	  if (NEWPID()->prolongation)
	    prolo = true;
	if (prolo) {
	  /* if we are too refined, we unrefine before sending */
	  cell.flags |= leaf;
	  array_append (a, &cell, sizeof(Cell));
	  cell.flags &= ~leaf;
	}
	else
	  array_append (a, &cell, sizeof(Cell));
      }
      else
	array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  }

  if (empty)
    a->len = 0;
  return a;
}

@def foreach_tree(t, size, list)
{
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);
  scalar * _list = list;
  char * _i = (char *) (t)->p;
  foreach_cell_all() {
    Cell * c = (Cell *) _i;
    if (c->flags & _sent) {
      _i += size;
@

@def end_foreach_tree()
    }
    else
      _i += sizeof(Cell);
    if (c->flags & _next) {
      assert (c->neighbors);
      if (!(c->flags & leaf) && is_leaf(cell) &&
	  (!is_newpided() || !NEWPID()->leaf))
	/* refined */
	refine_cell (point, _list, 0, NULL);
      else if (!cell.neighbors)
	/* prolongation */
	alloc_children (point);
    }
    else
      continue;
  } end_foreach_cell_all();
}
@

Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
  foreach_cell() {
    // root cells
    bool root = false;
    if ((!is_local(cell) || NEWPID()->pid - 1 != nextpid) && is_refined(cell)) {
      foreach_child()
	if (is_local(cell) && NEWPID()->pid - 1 == nextpid)
	  root = true, break;
      if (root && cell.pid != nextpid) {
	foreach_neighbor()
	  if (cell.pid != nextpid && is_newpided()) {
	    if (fp)
	      fprintf (fp, "%g %g %g %d %d root\n",
		       x, y, z, NEWPID()->pid - 1, cell.pid);
	    cell.flags |= sent;
	  }
      }
    }
    // children
    if ((is_local(cell) && NEWPID()->pid - 1 == nextpid) || root) {
      foreach_neighbor(1)
	if (cell.neighbors && cell.pid != nextpid)
	  foreach_child()
	    if (cell.pid != nextpid && is_newpided()) {
	      if (fp)
		fprintf (fp, "%g %g %g %d %d nextpid\n",
			 x, y, z, NEWPID()->pid - 1, cell.pid);
	      cell.flags |= sent;
	    }
    }
    if (is_leaf(cell))
      continue;
  }

  return tree (sizeof(Cell) + datasize, newpid);
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

static bool receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, MOVED_TAG(),
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = malloc (a.len);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, MOVED_TAG(),
		    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");
    //    const unsigned short next = 1 << (user + 1);
    foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
	      datasize);
      assert (NEWPID()->pid > 0);
      if (fp)
	fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
		 x, y, z, NEWPID()->pid - 1, cell.pid,
		 c->flags & leaf,
		 cell.flags & leaf, from, NEWPID()->leaf);
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

trace
bool balance (double imbalance)
{
  assert (sizeof(NewPid) == sizeof(double));

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
  
  scalar newpid[];
  z_indexing (newpid, true);

  FILE * fp = NULL;
#ifdef DEBUGCOND
  extern double t;
  if (DEBUGCOND)
    fp = lfopen ("bal", "w");
#elif DEBUG_MPI
  fp = lfopen ("bal", "w");
#endif

  // compute new pid, stored in newpid[]
  foreach_cell_all() {
    if (is_local(cell)) {
      int pid = balanced_pid (newpid[], nt);
      NEWPID()->pid = clamp (pid, cell.pid - 1, cell.pid + 1) + 1;
      NEWPID()->leaf = is_leaf(cell);
      NEWPID()->prolongation = is_prolongation(cell);
      if (fp)
	fprintf (fp, "%g %g %d %d newpid\n", x, y, NEWPID()->pid - 1, cell.pid);
    }
    else
      newpid[] = 0;
  }
  for (int l = 0; l <= depth(); l++)
    boundary_iterate (restriction, {newpid}, l);

#ifdef DEBUGCOND
  extern double t;
  if (DEBUGCOND) {
    char name[80];
    sprintf (name, "colls-before-%d", pid());
    FILE * fp = fopen (name, "w");
    output_cells (fp);
    fclose (fp);

    sprintf (name, "pid-before-%d", pid());
    fp = fopen (name, "w");
    foreach_cell() {
      fprintf (fp, "%g %g %g %d %d %d %d\n",
	       x, y, z, cell.pid, NEWPID()->pid - 1,
	       NEWPID()->leaf, is_leaf(cell));
      if (NEWPID()->leaf)
	assert (is_leaf(cell));
    }
    fclose (fp);  
  }
#endif // DEBUGCOND
  
  Array * anext = neighborhood (newpid, pid() + 1, fp);
  Array * aprev = neighborhood (newpid, pid() - 1, NULL);

  if (fp)
    fflush (fp);
  
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
  if (pid() < npe() - 1 && receive_tree (pid() + 1, newpid, fp))
    quadtree->dirty = true;
  if (pid() > 0 &&         receive_tree (pid() - 1, newpid, fp))
    quadtree->dirty = true;

  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  foreach_cell() {
    if (is_newpided() && NEWPID()->leaf && !is_leaf(cell)) {
      assert (cell.pid >= 0);
      assert (cell.neighbors);
      coarsen_cell_recursive (point, NULL);
    }
    if (is_leaf(cell))
      continue;
  }
  
  foreach_cell_all() {
    if (fp)
      fprintf (fp, "%g %g %g %d %d %d %d nowpid\n",
	       x, y, z, NEWPID()->pid - 1, cell.pid, is_leaf(cell),
	       NEWPID()->leaf);
    if (level < depth() && !cell.neighbors &&
	point.i > 0 && point.i <= (1 << level) + 2 &&
#if dimension >= 2
	point.j > 0 && point.j <= (1 << level) + 2 &&
#endif
#if dimension >= 3
	point.k > 0 && point.k <= (1 << level) + 2 &&
#endif	
	allocated_child(0)) {
      if (fp)
	fprintf (fp, "%g %g %g %d %d freechildren\n",
		 x, y, z, NEWPID()->pid - 1, cell.pid);
      free_children (point);
    }
  }

  if (fp)
    fflush (fp);
  
  // set new pids
  int pid_changed = false;
  foreach_cell() {
    if (is_newpided()) {
      if (fp)
	fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
		 x, y, z, NEWPID()->pid - 1, cell.pid,
		 is_leaf(cell), cell.neighbors, NEWPID()->leaf);
      if (cell.pid != NEWPID()->pid - 1) {
	cell.pid = NEWPID()->pid - 1;
	cell.flags &= ~(active|border);
	if (is_local(cell))
	  cell.flags |= active;
	pid_changed = true;
      }
    }
    else if (level > 0 && ((NewPid *)&coarse(newpid,0))->leaf)
      cell.pid = aparent(0).pid;
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
  
  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update();
  
  return pid_changed;
}
