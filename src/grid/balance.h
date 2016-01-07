// Dynamic load-balancing

void balance()
{
  int nf = 0;
  foreach(reduction(+:nf))
    nf++;
  nf = max(1, nf/npe());

  // fixme: check balancing
  // fixme: check that pid() - 1, pid() + 1 are accessible
  
  scalar index[];
  z_indexing (index, true);
  static const int moved = 1 << user;
  foreach() {
    int pid = min(npe() - 1, ((int)index[])/nf);
    index[] = clamp(pid, pid() - 1, pid() + 1);
    if (index[] != pid())
      foreach_neighbor()
	if (!is_boundary(cell)) {
	  cell.flags |= moved;
	  if (level > 0)
	    aparent(0).flags |= moved;
	}
  }
  
  // restriction
  boundary_iterate (restriction, {index}, depth());
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level(l) {
      index[] = fine(index,0,1,1);
      if (is_local(cell) && index[] != pid()) {
	foreach_neighbor()
	  if (!is_boundary(cell)) {
	    cell.flags |= moved;
	    if (level > 0)
	      aparent(0).flags |= moved;
	  }
      }
      else if (level > 0 && (cell.flags & moved))
	aparent(0).flags |= moved;
    }
    boundary_iterate (restriction, {index}, l);
  }

  // linear quadtree
  Array * a = array_new (sizeof(Cell));
  static const int stop = 1 << (user + 1);
  foreach_cell() {
    if (cell.flags & moved) {
      if (!is_local(cell))
	cell.flags &= ~moved;
      else if (index[] != pid()) {
	quadtree->dirty = true;
	cell.flags &= ~(active|border);
	cell.pid = index[];
	if (is_leaf(cell) && cell.neighbors > 0) {
	  int pid = cell.pid;
	  foreach_child()
	    cell.pid = pid;
	}
      }
      array_append (a, &cell);
      cell.flags &= ~moved;
    }
    else { // not moved
      short flags = cell.flags;
      if (is_local(cell))
	cell.flags |= moved;
      cell.flags |= stop;
      array_append (a, &cell);
      cell.flags = flags;
      continue;
    }
    if (is_leaf(cell))
      continue;
  }

  /* Send halo mesh for each neighboring process. */
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  RcvPid * snd = mpi->restriction.snd;
  MPI_Request r[2*snd->npid];
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid;
    MPI_Isend (&a->len, 1, MPI_INT, pid, MOVED_TAG(), MPI_COMM_WORLD, &r[2*i]);
    MPI_Isend (a->p, a->len*sizeof(Cell),
	       MPI_BYTE, pid, MOVED_TAG(), MPI_COMM_WORLD, &r[2*i+1]);
  }

  /* Receive halo mesh from each neighboring process. */
  RcvPid * rcv = mpi->restriction.rcv;
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid;
    // we could use MPI_Probe, but apparently it's not a good idea,
    // see http://cw.squyres.com/ and
    // http://cw.squyres.com/columns/2004-07-CW-MPI-Mechanic.pdf
    int len;
    MPI_Recv (&len, 1, MPI_INT, pid, MOVED_TAG(), MPI_COMM_WORLD,
	      MPI_STATUS_IGNORE);
    Cell a[len], * c = a;
    MPI_Recv (a, len*sizeof(Cell), MPI_BYTE, pid, MOVED_TAG(), MPI_COMM_WORLD,
	      MPI_STATUS_IGNORE);
    foreach_cell() {
      if (!(c->flags & leaf) && is_leaf(cell))
	refine_cell (point, NULL, 0, NULL);
      if ((c->flags & moved) && cell.pid != c->pid) {
	quadtree->dirty = true;
	cell.pid = c->pid;
	if (is_leaf(cell)) {
	  if (cell.neighbors > 0) {
	    int pid = cell.pid;
	    foreach_child()
	      cell.pid = pid;
	  }
	  if (is_local(cell))
	    cell.flags |= active;
	  else // not local
	    cell.flags &= ~(active|border);
	}
	else if (is_local(cell))
	  cell.flags |= active;
      }
      if ((c++->flags & stop) || is_leaf(cell))
	continue;
    }
  }
  
  /* check that ghost values were received OK and free send buffers */
  MPI_Waitall (2*snd->npid, r, MPI_STATUSES_IGNORE);
  array_free (a);

  // update active cells
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
  
  flag_border_cells();
  
  mpi_boundary_update();  
}
