// Dynamic load-balancing

static inline void update_pid (Point point, int pid)
{
  if (cell.neighbors > 0 && is_leaf(cell))
    foreach_child()
      cell.pid = pid;
  cell.pid = pid;
  cell.flags &= ~(active|border|ignored);
  if (is_local(cell))
    cell.flags |= active;
}

static bool send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_INT, to, MOVED_TAG(), MPI_COMM_WORLD, &r[0]);
  if (a->len > 1) {
    MPI_Isend (a->p, a->len*a->size, MPI_BYTE, to, MOVED_TAG(),
	       MPI_COMM_WORLD, &r[1]);
    return true;
  }
  else
    return false;
}

static bool receive_tree (size_t size, int from, FILE * fp)
{
  int len;
  MPI_Recv (&len, 1, MPI_INT, from, MOVED_TAG(),
	    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (len > 1) {
    int flag = 1 << user, pflag = 1 << (user + 1);
    char a[len*size], * i = a;
    if (fp)
      fprintf (fp, "receiving %d from %d\n", len, from);
    MPI_Recv (a, len*size, MPI_BYTE, from, MOVED_TAG(),
	      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    foreach_cell() {
      Cell * c = (Cell *) i; i += size;
      if (c->flags & flag) {
	if (fp)
	  fprintf (fp, "recv %g %g %d %d %d %d\n",
		   x, y, c->pid, cell.pid, from, c->flags & active);
	update_pid (point, c->pid);
	if ((c->flags & active) && is_leaf(cell) && is_local(cell))
	  /* fixme: this could be optimised by only sending data for
	     local leaf cells */
	  memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
		  datasize);
      }
      if (c->flags & leaf)
	continue;
      else if (is_leaf(cell))
	refine_cell (point, NULL, 0, NULL);
      if (!(c->flags & pflag))
	continue;
    }
    return true;
  }
  else
    return false;
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 1)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

bool indexing (scalar index, double imbalance)
{
  long nl = 0;
  foreach()
    nl++;
  long nt = nl, nmin = nl, nmax = nl;
  // fixme: do all reductions in one go
  mpi_all_reduce (nt,   MPI_LONG, MPI_SUM);
  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  long ne = max(1, nt/npe());
  if (nmax - nmin <= nt - ne*npe() || nmax - nmin <= ne*imbalance)
    return false;
  
  //  quadtree_trash ({index});
  z_indexing (index, true);
  // parallel index restriction
  for (int l = depth(); l >= 0; l--)
    mpi_boundary_restriction_children (mpi_boundary, {index}, l);

  foreach_cell() {
    if (cell.flags & ignored)
      index[] = undefined;
    else {
#if 1     
      int pid = min(npe() - 1, ((int)index[])/ne);
      index[] = clamp (pid, cell.pid - 1, cell.pid + 1);
#else
      index[] = cell.pid;
#endif
    }
    if (is_leaf(cell))
      continue;
  }

  return true;
}

/**
   Returns a linear quadtree which contains cells whose 'index[]' is
   equal to 'pid' and their neighborhood. */

Array * neighborhood (scalar index, int pid, size_t size)
{
  static const int flag = 1 << user, pflag = 1 << (user + 1);
  foreach_cell_post (!is_leaf(cell)) {
    //    assert (!(cell.flags & (flag|pflag)));
    if (is_local(cell) && index[] == pid) {
      //      fprintf (fp, "%g %g\n", x, y);
      cell.flags |= (flag|pflag);
    }
    else if (cell.pid != pid) {
      // =========================================
      // non-local cell: do we need to receive it?
      bool used = false;
      foreach_neighbor()
	if (allocated(0) && !is_boundary(cell)) {
	  if (!(cell.flags & ignored) && !is_prolongation(cell) &&
	      is_local(cell) && index[] == pid)
	    used = true, break;
	  // root cell: fixme: optimise?
	  if (is_refined(cell))
	    foreach_child()
	      if (is_local(cell) && index[] == pid)
		used = true, break;
	  if (used)
	    break;
	  /* see the corresponding condition in
	     mpi_boundary_refine(), and src/test/mpi-refine1.c
	     on why `level > 0 && is_local(aparent(0)))` is necessary. */
	  if (level > 0 && is_local(aparent(0)) && coarse(index,0) == pid)
	    used = true, break;
	}
      if (used) {
	cell.flags |= (flag|pflag);
	if (!is_leaf(cell))
	  foreach_child()
	    cell.flags |= (flag|pflag);
      }
    }
    if (level > 0 && (cell.flags & pflag))
      aparent(0).flags |= pflag;
  }

  Array * a = array_new (size);
  foreach_cell() {
    bool stop = !(cell.flags & pflag);
    if (cell.flags & flag) {
      int pid = cell.pid; cell.pid = index[];
      array_append (a, &cell);
      cell.pid = pid;
    }
    else
      array_append (a, &cell);
    cell.flags &= ~(flag|pflag);
    if (stop || is_leaf(cell))
      continue;
  }
  
  return a;
}
 
bool balance (double imbalance)
{
  scalar index[];
  if (!indexing (index, imbalance))
    return false;

  char name[80];
  sprintf (name, "bal-%d", pid());
  FILE * fp = fopen (name, "w");

#if 1  
  foreach_cell() {
    fprintf (fp, "%g %g %g %d index\n", x, y, index[], !(cell.flags & ignored));
    if (is_leaf(cell))
      continue;
  }
  fflush (fp);
#endif

  debug_mpi (NULL);
  
  Array * aprev = neighborhood (index, pid() - 1, sizeof(Cell) + datasize);
  Array * anext = neighborhood (index, pid() + 1, sizeof(Cell) + datasize);
  
  /* fixme: this should be optimisable using the fact that "send to
     next" must be mutually exclusive with "receive from next" */
  // send mesh to previous/next process
  MPI_Request rprev[2], rnext[2];
  if (pid() > 0         && send_tree (aprev, pid() - 1, rprev))
    quadtree->dirty = true;
  if (pid() < npe() - 1 && send_tree (anext, pid() + 1, rnext))
    quadtree->dirty = true;

  // update pids
  foreach_cell() {
    if (!(cell.flags & ignored) && cell.pid != index[])
      update_pid (point, index[]);
    if (is_leaf(cell))
      continue;
  }
  
  // receive mesh from next/previous process
  if (pid() < npe() - 1 && receive_tree (anext->size, pid() + 1, fp))
    quadtree->dirty = true;
  if (pid() > 0 &&         receive_tree (aprev->size, pid() - 1, fp))
    quadtree->dirty = true;
  
  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);
  
  if (1/*quadtree->dirty*/) {
#if 1
    // update active cells: can this be done above
    foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
	bool inactive = true;
	foreach_child()
	  if (is_active(cell))
	    inactive = false, break;
	if (inactive)
	  cell.flags &= ~active;
	else
	  cell.flags |= active;
      }
#endif  
    flag_border_cells (fp); // can this be done above?
  }
    
  //  fprintf (stderr, "pid: %d changed: %d\n", pid(), changed);
    
  fclose (fp);

  MPI_Barrier (MPI_COMM_WORLD);
  
  if (true/*pid_changed*/)
    mpi_boundary_update();

  return true;
}
