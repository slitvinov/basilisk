#define DEBUG 0

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
  double * buf;      // MPI buffer
  MPI_Request r;     // MPI request 
} Rcv;

typedef struct {
  int nrcv;   // the number of PEs to receive ghost values from
  Rcv * rcv;  // ghost cells to receive

  int nsnd;   // the number of PEs to send ghost values to
  Rcv * snd;  // ghost cells to send

  bool children;   // whether to send children or only parent
} SndRcv;

typedef struct {
  Boundary parent;
  
  SndRcv restriction;
  SndRcv prolongation;

  int * adjacency; // adjacency[x] is different from zero if we need to send
                   // data to proc x
} MpiBoundary;

static void rcv_append (Point point, Rcv * rcv)
{
  if (level > rcv->depth) {
    rcv->halo = realloc (rcv->halo, (level + 1)*sizeof (CacheLevel));
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
}

void rcv_print (Rcv * rcv)
{
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      foreach_cache_level(rcv->halo[i], i,)
	fprintf (stderr, "%g %g %d %d\n", x, y, level, rcv->pid);
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    free (rcv->buf);
    rcv->buf = NULL;
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      free (rcv->halo[i].p);
  free (rcv->halo);
}

static void snd_rcv_append (SndRcv * m, Point point, int pid)
{
  int i;
  for (i = 0; i < m->nrcv; i++)
    if (pid == m->rcv[i].pid)
      break;
  if (i == m->nrcv) {
    m->rcv = realloc (m->rcv, ++m->nrcv*sizeof (Rcv));
    Rcv * rcv = &m->rcv[m->nrcv-1];
    rcv->pid = pid;
    rcv->depth = 0;
    rcv->halo = malloc (sizeof (CacheLevel));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  if (level > 0)
    aparent(0,0).pid = pid;
  rcv_append (point, &m->rcv[i]);
}

static void snd_rcv_sync (SndRcv * m, scalar * list, int l)
{
  prof_start();
  
  // fixme: 4 should be 2**dimension
  int len = list_len (list)*(m->children ? 4 : 1); 

  /* send ghost values */
  for (int i = 0; i < m->nsnd; i++) {
    Rcv * snd = &m->snd[i];
    if (l <= snd->depth && snd->halo[l].n > 0) {
      rcv_free_buf (snd);
      snd->buf = malloc (sizeof (double)*snd->halo[l].n*len);
      double * b = snd->buf;
      if (m->children)
	foreach_cache_level(snd->halo[l], l,)
	  foreach_child() {
	    assert (allocated(0,0));
	    for (scalar s in list)
	      *b++ = s[];
	  }
      else
	foreach_cache_level(snd->halo[l], l,)
	  for (scalar s in list)
	    *b++ = s[];
      MPI_Isend (snd->buf, snd->halo[l].n*len, MPI_DOUBLE, snd->pid, l, 
		 MPI_COMM_WORLD, &snd->r);
    }
  }
    
  /* receive ghost values */
  for (int i = 0; i < m->nrcv; i++) {
    Rcv rcv = m->rcv[i];
    if (l <= rcv.depth && rcv.halo[l].n > 0) {
      double buf[rcv.halo[l].n*len], * b = buf;
      MPI_Status s;
      MPI_Recv (buf, rcv.halo[l].n*len, MPI_DOUBLE, rcv.pid, l, 
		MPI_COMM_WORLD, &s);
      if (m->children)
	foreach_cache_level(rcv.halo[l], l,)
	  foreach_child()
	    for (scalar s in list)
	      s[] = *b++;
      else
	foreach_cache_level(rcv.halo[l], l,)
	  for (scalar s in list)
	    s[] = *b++;
    }
  }
    
  prof_stop();
}

void snd_rcv_print (SndRcv * m, FILE * fp)
{
  for (int i = 0; i < m->nrcv; i++)
    for (int j = 0; j <= m->rcv[i].depth; j++)
      foreach_cache_level (m->rcv[i].halo[j], j,)
	fprintf (fp, "%g %g %d %d %d\n",
		 x, y, m->rcv[i].pid, cell.neighbors, level);
}

void rcv_graph (SndRcv * m, FILE * fp)
{
  for (int i = 0; i < m->nrcv; i++) {
    fprintf (fp, "%d -> %d [label=\"", m->rcv[i].pid, pid());
    for (int j = 0; j <= m->rcv[i].depth; j++) {
      int n = 0;
      foreach_cache_level (m->rcv[i].halo[j], j,)
	n++;
      fprintf (fp, "%d ", n);
    }
    fprintf (fp, "\"];\n");
  }
}

static void snd_free (SndRcv * m)
{
  for (int i = 0; i < m->nsnd; i++)
    rcv_destroy (&m->snd[i]);
  free (m->snd);
  m->snd = NULL; m->nsnd = 0;
}

static void rcv_free (SndRcv * m)
{
  for (int i = 0; i < m->nrcv; i++)
    rcv_destroy (&m->rcv[i]);
  free (m->rcv);
  m->rcv = NULL; m->nrcv = 0;
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_free (m);
  snd_free (m);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->restriction);
  snd_rcv_destroy (&m->prolongation);
  free (m->adjacency);
  free (m);
}

static void mpi_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_sync (&m->restriction, list, l);
}

static void mpi_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, int l, int depth)
{
  if (l > 0) {
    MpiBoundary * m = (MpiBoundary *) b;
    snd_rcv_sync (&m->prolongation, list, l - 1);
  }
}

static Boundary * mpi_boundary = NULL;

void mpi_boundary_new()
{
  mpi_boundary = calloc (1, sizeof (MpiBoundary));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->restriction = mpi_boundary->halo_restriction = 
    mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation =
    mpi_boundary_halo_prolongation;
  MpiBoundary * b = (MpiBoundary *) mpi_boundary;
  b->adjacency = calloc (npe(), sizeof (int));
  b->prolongation.children = true;
  add_boundary (mpi_boundary);
}

static void cache_level_send (CacheLevel * halo, int dest, int tag)
{
  MPI_Request r;
  int n = halo ? halo->n : 0;
  MPI_Isend (&n, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD, &r);
  if (n > 0)
    MPI_Isend (halo->p, 2*halo->n, MPI_INT, dest, tag, MPI_COMM_WORLD, &r);
}

static void cache_level_receive (CacheLevel * halo, int src, int tag)
{
  MPI_Status s;
  MPI_Recv (&halo->n, 1, MPI_INT, src, tag + 1, MPI_COMM_WORLD, &s);
  if (halo->n > 0) {
    halo->p = malloc (halo->n*sizeof(IndexLevel));
    MPI_Recv (halo->p, 2*halo->n, MPI_INT, src, tag, MPI_COMM_WORLD, &s);
  }
}
	
void debug_mpi()
{
  void output_cells (FILE * fp);

  char name[80];
  sprintf (name, "halo-%d", pid());
  FILE * fp = fopen (name, "w");

  // local halo
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l)
      foreach_child()
        fprintf (fp, "%g %g %d\n", x, y, level);
  fclose (fp);

  sprintf (name, "cells-%d", pid());
  fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  sprintf (name, "neighbors-%d", pid());
  fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %d\n", x, y, cell.neighbors);
  fclose (fp);

  // local restriction
  sprintf (name, "restriction-%d", pid());
  fp = fopen (name, "w");
  for (int l = 0; l < depth(); l++)
    foreach_halo (restriction, l)
      fprintf (fp, "%g %g %d %d\n", x, y, level, cell.neighbors);
  fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  
  sprintf (name, "mpi-restriction-%d", pid());
  fp = fopen (name, "w");
  snd_rcv_print (&mpi->restriction, fp);
  fclose (fp);

  sprintf (name, "mpi-prolongation-%d", pid());
  fp = fopen (name, "w");
  snd_rcv_print (&mpi->prolongation, fp);
  fclose (fp);
}

static void mpi_boundary_sync (MpiBoundary * mpi)
{
  prof_start();

  SndRcv * restriction = &mpi->restriction;
  SndRcv * prolongation = &mpi->prolongation;

  /* we send the ghost cell indices to their respective PEs */
  for (int i = 0; i < restriction->nrcv; i++) {
    Rcv res = restriction->rcv[i], * pro = NULL;
    /* pro is the prolongation received from the same pid */
    for (int j = 0; j < prolongation->nrcv && !pro; j++)
      if (prolongation->rcv[j].pid == res.pid)
	pro = &prolongation->rcv[j];

    MPI_Request r;
    MPI_Isend (&res.depth, 1, MPI_INT, res.pid, 0, MPI_COMM_WORLD, &r);
    for (int j = 0; j <= res.depth; j++) {
      // restriction buffer
      cache_level_send (&res.halo[j], res.pid, j);
      // prolongation buffer
      if (pro && j <= pro->depth)
	cache_level_send (&pro->halo[j], res.pid, j);
      else
	cache_level_send (NULL, res.pid, j);
    }
  }

  /* we free the old buffers */
  snd_free (restriction);
  snd_free (prolongation);

  /* we receive the ghost cell indices from their respective PEs */
  for (int i = 0; i < npe(); i++)
    if (mpi->adjacency[i]) {
      restriction->snd = realloc (restriction->snd,
				  ++restriction->nsnd*sizeof (Rcv));
      prolongation->snd = realloc (prolongation->snd,
				  ++prolongation->nsnd*sizeof (Rcv));
      Rcv * res = &restriction->snd[restriction->nsnd-1];
      Rcv * pro = &prolongation->snd[prolongation->nsnd-1];
      res->pid = pro->pid = i;
      MPI_Status s;
      MPI_Recv (&res->depth, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &s);
      pro->depth = res->depth;
      pro->buf = res->buf = NULL;
      res->halo = malloc ((res->depth + 1)*sizeof (CacheLevel));
      pro->halo = malloc ((pro->depth + 1)*sizeof (CacheLevel));
      for (int j = 0; j <= res->depth; j++) {
	cache_level_receive (&res->halo[j], res->pid, j);
	cache_level_receive (&pro->halo[j], res->pid, j);
      }
    }
  
  prof_stop();
  
#if DEBUG
  debug_mpi();
#endif  
}

void mpi_boundary_update (MpiBoundary * m)
{
  prof_start();

  rcv_free (&m->restriction);
  rcv_free (&m->prolongation);

  foreach_cell() {
    if (!is_active(cell)) {
      bool neighbors = false;
      /* If a leaf has (active) refined neighbors, we add it to the
	 prolongation buffer. */
      if (is_leaf(cell))
	for (int k = -GHOSTS/2; k <= GHOSTS/2 && !neighbors; k++)
	  for (int l = -GHOSTS/2; l <= GHOSTS/2 && !neighbors; l++)
	    if (allocated(k,l) && is_active(neighbor(k,l)) &&
		!is_leaf(neighbor(k,l))) {
	      snd_rcv_append (&m->prolongation, point, cell.pid);
	      neighbors = true;
	    }
      /* If a cell is in the neighborhood of an active cell, we add it
	 to the restriction buffer. */
      for (int k = -GHOSTS; k <= GHOSTS && !neighbors; k++)
	for (int l = -GHOSTS; l <= GHOSTS && !neighbors; l++)
	  if (allocated(k,l) && is_active(neighbor(k,l)))
	    neighbors = true;
      if (neighbors) {
	snd_rcv_append (&m->restriction, point, cell.pid);
	cell.flags |= halo;
      }
      else
	cell.flags &= ~halo;
    }
    if (is_leaf(cell))
      continue;
  }

  prof_stop();
}
  
/*
  Returns a linear quadtree containing the leaves which are in the
  (5x5) neighborhood of cells belonging to the (remote) process 'pid'.
 */
Array * remote_leaves (unsigned pid)
{
  static const int remote = 1 << user;
  foreach_cell_post (is_active (cell))
    if (is_active(cell)) {
      if (is_leaf(cell)) {
	bool neighbors = false;
	if (is_leaf(cell))
	  for (int k = -GHOSTS; k <= GHOSTS && !neighbors; k++)
	    for (int l = -GHOSTS; l <= GHOSTS && !neighbors; l++)
	      if (neighbor(k,l).pid == pid) {
		cell.flags |= remote;
		neighbors = true;
	      }
      }
      if ((cell.flags & remote) && level > 0)
	aparent(0,0).flags |= remote;
    }

  Array * a = array_new (sizeof(unsigned));
  foreach_cell() {
    unsigned flags = cell.flags;
    if (flags & remote) {
      flags &= ~remote;
      cell.flags = flags;
    }
    else
      flags |= leaf;
    array_append (a, &flags);
    if (flags & leaf)
      continue;
  }
  return a;
}

/**
   Match the refinement given by the linear quadtree 'a'. 
   Returns the number of active cells refined.
*/
int match_refine (unsigned * a, scalar * list)
{
  int nactive = 0;
  if (a != NULL)
    foreach_cell() {
      unsigned flags = *a++;
      if (flags & leaf)
	continue;
      else if (is_leaf(cell))
	point = refine_cell (point, list, 0, &nactive);	
    }
  return nactive;
}

void mpi_boundary_match (MpiBoundary * mpi)
{
  prof_start();
  
  /* Send halo mesh for each neighboring process. */
  SndRcv * m = &mpi->restriction;
  Array * a[m->nsnd];
  MPI_Request r[m->nsnd];
  for (int i = 0; i < m->nsnd; i++) {
    Rcv snd = m->snd[i];
    a[i] = remote_leaves (snd.pid);
    //    fprintf (stderr, "  remote %d %ld\n", i, a[i]->len);
    assert (a[i]->len > 0);
    MPI_Isend (&a[i]->len, 1, MPI_INT, snd.pid, 0, MPI_COMM_WORLD, &r[i]);
    MPI_Isend (a[i]->p, a[i]->len, MPI_INT, snd.pid, 0, MPI_COMM_WORLD, &r[i]);
  }
    
  /* Receive halo mesh from each neighboring process. */
  int refined = 0;
  for (int i = 0; i < m->nrcv; i++) {
    Rcv rcv = m->rcv[i];
    MPI_Status s;
    int len;
    MPI_Recv (&len, 1, MPI_INT, rcv.pid, 0, MPI_COMM_WORLD, &s);
    unsigned * a = malloc (len*sizeof(unsigned));
    MPI_Recv (a, len, MPI_INT, rcv.pid, 0, MPI_COMM_WORLD, &s);
    refined += match_refine (a, NULL);
    free (a);
  }
  
  /* check that ghost values were received OK and free send buffers */
  for (int i = 0; i < m->nsnd; i++) {
    MPI_Status s;
    MPI_Wait (&r[i], &s);
    array_free (a[i]);
  }
  
  prof_stop();
  
  /* If an active cell has been refined (due to the 2:1 refinement
     constraint), we need to repeat the process. */
  mpi_all_reduce (refined, MPI_INT, MPI_SUM);
  if (refined)
    mpi_boundary_match (mpi);
}

void mpi_boundary_refine (void * p, scalar * list)
{
  MpiBoundary * mpi = p ? p : (MpiBoundary *) mpi_boundary;
  mpi_boundary_match (mpi);
  mpi_boundary_update (mpi);
  mpi_boundary_sync (mpi);
}

void mpi_partitioning()
{
  prof_start();

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  int nf = 0;
  foreach()
    nf++;
  nf = max(1, nf/npe());

  /* set the pid of each cell */
  int i = 0;
  foreach_cell_post (is_active (cell) || cell.neighbors > 0)
    if (is_active (cell)) {
      if (is_leaf (cell)) {
	cell.pid = min(npe() - 1, i/nf); i++;
	if (cell.neighbors > 0)
	  for (int i = 0; i <= 1; i++)
	    for (int j = 0; j <= 1; j++)
	      child(i,j).pid = cell.pid;
      }
      else {
	cell.pid = -1;
	for (int i = 0; i <= 1; i++)
	  for (int j = 0; j <= 1; j++)
	    if (cell.pid == -1)
	      cell.pid = child(i,j).pid;
	    else if (child(i,j).pid == pid())
	      cell.pid = pid();
      }
      if (cell.pid != pid())
	cell.flags &= ~active;
    }

  /* adjacency[i] is different from zero if we need to receive data
     from proc i */
  int adjacency[npe()];
  for (int i = 0; i < npe(); i++)
    adjacency[i] = 0;
  
  /* we build arrays of ghost cell indices for restriction and prolongation */
  ((Quadtree *)grid)->dirty = true;
  foreach_cell_post (!is_leaf(cell))
    if (!is_active(cell)) {
      int neighbors = 0;
      for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	  if (allocated(k,l)) {
	    /* If the cell is a leaf, we check whether it has active
	       and refined neighbors. */
	    if (is_leaf(cell)) {
	      if (is_active(neighbor(k,l)) && !is_leaf(neighbor(k,l)))
		neighbors++;
	    }
	    /* We update the number of neighbors only for non-active
	       (i.e. remote) cells. */
	    else if (!is_active(neighbor(k,l))) {
	      neighbor(k,l).neighbors--;
	      if (neighbor(k,l).neighbors == 0)
		free_children (point,k,l);
	    }
	  }

      /* If the cell is a leaf and has active and refined neighbors,
	 we add it to the prolongation buffer. */
      if (neighbors)
	snd_rcv_append (&m->prolongation, point, cell.pid);
      else
	/* If the cell does not have active and refined neighbors, 
	   is it in the neighborhood of an active cell? */
	for (int k = -GHOSTS; k <= GHOSTS && !neighbors; k++)
	  for (int l = -GHOSTS; l <= GHOSTS && !neighbors; l++)
	    if (allocated(k,l) && is_active(neighbor(k,l)))
	      neighbors = true;
      
      /* If the cell has active and refined neighbors, or it is in the
	 neighborhood of an active cell, then it is a restriction. */
      if (neighbors) {
	snd_rcv_append (&m->restriction, point, cell.pid);
	adjacency[cell.pid]++;
      }
    }

  /* we do a global transpose of the adjacency vector
     i.e. m->adjacency[i] is different from zero if we need to send
     data to proc i */
  MPI_Alltoall (adjacency, 1, MPI_INT, m->adjacency, 1, MPI_INT,
		MPI_COMM_WORLD);
  
  prof_stop();
  
  mpi_boundary_sync (m);  
}

/**
# *z_indexing()*: fills *index* with the Z-ordering index.
   
On a single processor, we would just need something like

~~~literatec
double i = 0;
foreach()
  index[] = i++;
~~~

In parallel, this is a bit more difficult. */

void z_indexing (scalar index)
{
  /**
  ## Size of subtrees

  We first compute the size of each subtree. We use *index* to store
  this. */
  
  scalar size = index;

  /**
  The size of leaf "subtrees" is one. */

  foreach()
    size[] = 1;
  
  /**
  We do a (parallel) restriction to compute the size of non-leaf
  subtrees. */

  boundary_iterate (restriction, {size}, depth());
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level(l)
      size[] = (fine(size,0,0) + fine(size,1,0) + 
		fine(size,0,1) + fine(size,1,1));
    boundary_iterate (restriction, {size}, l);
  }

  /**
  ## Indexing
  
  Indexing can then be done locally. */

  double i = 0;
  foreach_cell() {
    if (is_leaf(cell)) {
      index[] = i++;
      continue;
    }
    else if (!is_active(cell)) {
      i += size[];
      continue;
    }
  }
}

