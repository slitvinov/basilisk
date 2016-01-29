// Dynamic load-balancing

static inline void update_pid (Point point, int pid)
{
  if (cell.pid == pid)
    return;
  if (is_leaf(cell) && cell.neighbors > 0)
    foreach_child()
      cell.pid = pid;
  cell.pid = pid;
  cell.flags &= ~(active|border|ignored);
  if (is_local(cell))
    cell.flags |= active;
}

void check_active()
{
  foreach_cell_post (!is_leaf (cell)) {
    if (is_local(cell))
      assert (is_active(cell));
    else if (!is_leaf(cell)) {
      short flags = cell.flags & ~active;
      foreach_child()
	if (is_active(cell))
	  flags |= active, break;
      if (flags & active)
	assert (is_active(cell));
      else
	assert (!is_active(cell));
    }
    else
      assert (!is_active(cell));
  }
}

static inline void mark_neighbor (Point point, FILE * fp, char t)
{
  static const short flag = 1 << user, pflag = 1 << (user + 1);
#if 0
  if (is_leaf(cell) && (level > 0 && cell.pid == aparent(0).pid)) {
    if (fp)
      fprintf (fp, "%g %g %c nux\n", x, y, t);
    aparent(0).flags |= flag;
    return;
  }
#endif
  if (fp)
    fprintf (fp, "%g %g %c neigh\n", x, y, t);
  cell.flags |= flag;
  if (!is_leaf(cell)) {
#if 0
    int pid = cell.pid;
    foreach_child()
      if (cell.pid != pid)
	pid = -1, break;
    if (pid < 0)
#endif
    {
      cell.flags |= pflag;
      foreach_child()
	cell.flags |= flag;
    }
  }
  if (level > 0)
    aparent(0).flags |= pflag;
}

/**
   Returns a linear quadtree which contains cells whose 'index[]' is
   equal to 'pid' and their neighborhood. */

Array * neighborhood (int pid, short locflag, FILE * fp)
{
  static const short flag = 1 << user, pflag = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
    foreach_cell() {
      if (level == l) {
	if (cell.flags & locflag) {
	  cell.flags |= flag;
	  foreach_neighbor()
	    if (!is_prolongation(cell) && !is_boundary(cell) &&
		!(cell.flags & locflag) && cell.pid != pid)
	      mark_neighbor (point, fp, 'a');
	  foreach_neighbor(1)
	    if (!is_boundary(cell) && is_refined(cell))
	      foreach_child()
		if (!is_boundary(cell) &&
		    !(cell.flags & locflag) && cell.pid != pid)
		  mark_neighbor (point, fp, 'b');
	}
	else if (!is_leaf(cell)) {
	  // is it a root cell?
	  bool root = false;
	  foreach_child()
	    if (cell.flags & locflag)
	      root = true, break;
	  if (root)
	    foreach_neighbor()
	      if (!is_prolongation(cell) && !is_boundary(cell) &&
		  !(cell.flags & locflag) && cell.pid != pid)
		mark_neighbor (point, fp, 'c');
	}
	if (level > 0 && (cell.flags & (pflag|flag)))
	  aparent(0).flags |= pflag;
	continue;
      }
      if (is_leaf(cell))
	continue;
    }

  Array * a = array_new();
  foreach_cell() {
    bool stop = !(cell.flags & pflag);
    if ((cell.flags & locflag) && is_leaf(cell))
      array_append (a, &cell, sizeof(Cell) + datasize);
    else
      array_append (a, &cell, sizeof(Cell));
    cell.flags &= ~(flag|pflag|locflag);
    if (stop)
      continue;
  }

  return a;
}

#define SIZE (sizeof(Cell) + datasize)

Array * neighborhood1 (int pid, short locflag, FILE * fp)
{
  static const short flag = 1 << user, pflag = 1 << (user + 1);
  foreach_cell_post (!is_leaf(cell)) {
    if (cell.flags & locflag)
      cell.flags |= flag;
    else if (cell.pid != pid) {
      // =========================================
      // non-local cell: do we need to receive it?
      bool used = false;
      double x1 = x, y1 = y;
      foreach_neighbor()
	if (allocated(0) && !is_boundary(cell)) {
	  if ((cell.flags & locflag) && !is_prolongation(cell)) {
	    if (fp)
	      fprintf (fp, "%g %g a neigh\n", x1, y1);
	    used = true, break;
	  }
	  // root cell: fixme: optimise?
	  if (is_refined(cell))
	    foreach_child()
	      if (cell.flags & locflag) {
		if (fp)
		  fprintf (fp, "%g %g b neigh\n", x1, y1);
		used = true, break;
	      }
	  if (used)
	    break;
	  /* see the corresponding condition in
	     mpi_boundary_refine(), and src/test/mpi-refine1.c
	     on why `level > 0 && is_local(aparent(0)))` is necessary. */
	  if (level > 0 && (aparent(0).flags & locflag)) {
	    if (fp)
	      fprintf (fp, "%g %g c neigh\n", x1, y1);
	    used = true, break;
	  }
	}
      if (used) {
	cell.flags |= flag;
	if (!is_leaf(cell)) {
	  cell.flags |= pflag;
	  foreach_child()
	    cell.flags |= flag;
	}
      }
    }
    if (level > 0 && (cell.flags & (pflag|flag)))
      aparent(0).flags |= pflag;
  }

  Array * a = array_new();
  foreach_cell() {
    bool stop = !(cell.flags & pflag);
    array_append (a, &cell, SIZE);
    cell.flags &= ~(flag|pflag|locflag);
    if (stop)
      continue;
  }
  
  return a;
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

static bool receive_tree (int from, short locflag, FILE * fp)
{
  long len;
  MPI_Recv (&len, 1, MPI_LONG, from, MOVED_TAG(),
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
		   x, y, c->pid, cell.pid, from, c->flags & locflag);
	update_pid (point, c->pid);
	if ((c->flags & locflag) && is_leaf(cell))
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

bool balance (double imbalance)
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
  //  if (ne < 1000) ne = 1000; // fixme: this does not work
  if (nmax - nmin <= nt - ne*npe() || nmax - nmin <= ne*imbalance)
    return false;
  
  scalar index[];
  z_indexing (index, true);
  // parallel index restriction
  for (int l = depth(); l >= 0; l--)
    mpi_boundary_restriction_children (mpi_boundary, {index}, l);

  static const short prev = 1 << (user + 2), next = 1 << (user + 3);
  long nprev = 0, nnext = 0;
  bool pid_changed = false;
  foreach_cell_post (!is_leaf(cell)) {
    if (!(cell.flags & ignored)) {
      int pid = min(npe() - 1, ((int)index[])/ne);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (cell.pid != pid) {
	if (is_local(cell)) {
	  if (pid == cell.pid - 1)
	    cell.flags |= prev, nprev++;
	  else
	    cell.flags |= next, nnext++;
	}
	update_pid (point, pid);
	pid_changed = true;
      }
    }
#if 0
    // update active cells
    if (!is_leaf(cell) && !is_local(cell)) {
      short flags = cell.flags & ~active;
      foreach_child()
	if (is_active(cell))
	  flags |= active, break;
      cell.flags = flags;
    }
#endif
  }

  //  check_active();
  
  char name[200];
  sprintf (name, "bal-%d", pid());
  FILE * fp = fopen (name, "w");

#if 1
  Array * aprev = nprev ?
    neighborhood (pid() - 1, prev, sizeof(Cell) + datasize, NULL) :
    array_new (sizeof(Cell) + datasize);

  timer t = timer_start();
  Array * anext = nnext ?
    neighborhood (pid() + 1, next, sizeof(Cell) + datasize, NULL) :
    array_new (sizeof(Cell) + datasize);
#if 0
  static int ni = 0;
  fprintf (stderr, "timings %d %g %ld %ld %ld %ld %ld\n",
	   ni++, timer_elapsed(t), anext->len, ne, nl, nnext, nprev);
  fflush (stderr);
#endif
#else
  Array * aprev = neighborhood1 (pid() - 1, prev, sizeof(Cell) + datasize, NULL);
  FILE * fp1;
#if 0
  sprintf (name, "neigh1-%d", pid());
  fp1 = fopen (name, "w");
  neighborhood1 (pid() + 1, next, sizeof(Cell) + datasize, fp1);
  fclose (fp1);
#endif
  
  sprintf (name, "neigh-%d", pid());
  fp1 = fopen (name, "w");
  Array * anext = neighborhood1 (pid() + 1, next, sizeof(Cell) + datasize, fp1);
  fclose (fp1);
  
#if 0
  sprintf (name,
	   "awk '{print $1,$2,%d}' neigh-%d | sort > res-%d;"
	   "awk '{print $1,$2,%d}' neigh1-%d | sort | uniq > nei-%d;"
	   "diff nei-%d res-%d > diff-%d",
	   pid(), pid(), pid(), pid(), pid(), pid(), pid(), pid(), pid());
  system (name);
  sprintf (name, "diff-%d", pid());
  fp1 = fopen(name, "r");
  fseek (fp1, 0L, SEEK_END);
  assert (!ftell(fp1));
  fclose (fp1);
#endif
#endif
  
  /* fixme: this should be optimisable using the fact that "send to
     next" must be mutually exclusive with "receive from next" */
  // send mesh to previous/next process
  MPI_Request rprev[2], rnext[2];
  if (pid() > 0         && send_tree (aprev, pid() - 1, rprev))
    quadtree->dirty = true;
  if (pid() < npe() - 1 && send_tree (anext, pid() + 1, rnext))
    quadtree->dirty = true;

  //  check_active();
  
  // receive mesh from next/previous process
  if (pid() < npe() - 1 && receive_tree (anext->size, pid() + 1, prev, fp))
    quadtree->dirty = true;
  if (pid() > 0 &&         receive_tree (aprev->size, pid() - 1, next, fp))
    quadtree->dirty = true;
  
  //  check_active();

  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);
  
  if (quadtree->dirty) {
#if 1
    // update active cells: can this be done above
    foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
	short flags = cell.flags & ~active;
	foreach_child()
	  if (is_active(cell))
	    flags |= active, break;
	cell.flags = flags;
      }
#endif  
    flag_border_cells (fp); // can this be done above?
  }

  MPI_Barrier (MPI_COMM_WORLD);
  fclose (fp);

  if (quadtree->dirty || pid_changed)
    mpi_boundary_update();
    
  return true;
}
