#define DEBUG 0

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
} Rcv;

typedef struct {
  Boundary parent;

  int nrcv;   // the number of PEs to receive ghost values from
  Rcv * rcv;  // ghost cells to receive

  int nsnd;   // the number of PEs to send ghost values to
  Rcv * snd;  // ghost cells to send

  int * adjacency; // adjacency[x] is different from zero if we need to send
                   // data to proc x
  int npe;         // the total number of PEs
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

static void rcv_destroy (Rcv * rcv)
{
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      free (rcv->halo[i].p);
  free (rcv->halo);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  for (int i = 0; i < m->nrcv; i++)
    rcv_destroy (&m->rcv[i]);
  free (m->rcv);
  for (int i = 0; i < m->nsnd; i++)
    rcv_destroy (&m->snd[i]);
  free (m->snd);
  free (m->adjacency);
  free (m);
}

static void mpi_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  int len = list_len (list);

  /* send ghost values */
  double * buf[m->nsnd];
  MPI_Request r[m->nsnd];
  for (int i = 0; i < m->nsnd; i++) {
    Rcv snd = m->snd[i];
    if (l <= snd.depth && snd.halo[l].n > 0) {
      buf[i] = malloc (sizeof (double)*snd.halo[l].n*len);
      double * b = buf[i];
      foreach_cache_level(snd.halo[l], l,)
	for (scalar s in list)
	  *b++ = s[];
      MPI_Isend (buf[i], snd.halo[l].n*len, MPI_DOUBLE, snd.pid, l, 
		 MPI_COMM_WORLD, &r[i]);
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
      foreach_cache_level(rcv.halo[l], l,) {
#if DEBUG
	fprintf (stderr, "%g %g rcv\n", x, y);
#endif
	for (scalar s in list)
	  s[] = *b++;
      }
    }
  }
    
  /* check that ghost values were received OK and free send buffers */
  for (int i = 0; i < m->nsnd; i++)
    if (l <= m->snd[i].depth && m->snd[i].halo[l].n > 0) {
      MPI_Status s;
      MPI_Wait (&r[i], &s);
      free (buf[i]);
    }
}

static Boundary * mpi_boundary = NULL;

void mpi_boundary_new()
{
  mpi_boundary = calloc (1, sizeof (MpiBoundary));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->restriction = mpi_boundary->halo_restriction = 
    mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation = none;
  MpiBoundary * b = (MpiBoundary *) mpi_boundary;
  MPI_Comm_size (MPI_COMM_WORLD, &b->npe);
  b->adjacency = calloc (b->npe, sizeof (int));
  add_boundary (mpi_boundary);
}

static void mpi_boundary_sync (MpiBoundary * m)
{
  /* we send the ghost cell indices to their respective PEs */
  for (int i = 0; i < m->nrcv; i++) {
    Rcv rcv = m->rcv[i];
    MPI_Request r;
    MPI_Isend (&rcv.depth, 1, MPI_INT, rcv.pid, 0, MPI_COMM_WORLD, &r);
    for (int j = 0; j <= rcv.depth; j++) {
      CacheLevel halo = rcv.halo[j];
      MPI_Isend (&halo.n, 1, MPI_INT, rcv.pid, j + 1, MPI_COMM_WORLD, &r);
      if (halo.n > 0)
	MPI_Isend (halo.p, 2*halo.n, MPI_INT, rcv.pid, j, MPI_COMM_WORLD, &r);
    }
  }

  /* we free the old buffers */
  for (int i = 0; i < m->nsnd; i++)
    rcv_destroy (&m->snd[i]);
  free (m->snd);
  m->snd = NULL; m->nsnd = 0;
  
  /* we receive the ghost cell indices from their respective PEs */
  for (int i = 0; i < m->npe; i++)
    if (m->adjacency[i]) {
      m->snd = realloc (m->snd, ++m->nsnd*sizeof (Rcv));
      Rcv * snd = &m->snd[m->nsnd-1];
      snd->pid = i;
      MPI_Status s;
      MPI_Recv (&snd->depth, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &s);
      snd->halo = malloc ((snd->depth + 1)*sizeof (CacheLevel));
      for (int j = 0; j <= snd->depth; j++) {
	CacheLevel * halo = &snd->halo[j];
	MPI_Recv (&halo->n, 1, MPI_INT, i, j + 1, MPI_COMM_WORLD, &s);
	if (halo->n > 0) {
	  halo->p = malloc (halo->n*sizeof(IndexLevel));
	  MPI_Recv (halo->p, 2*halo->n, MPI_INT, i, j, MPI_COMM_WORLD, &s);
	}
      }
    }

#if DEBUG
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
  sprintf (name, "restrict-%d", pid());
  fp = fopen (name, "w");
  for (int l = 0; l < depth(); l++)
    foreach_halo (restriction, l)
      fprintf (fp, "%g %g %d %d\n", x, y, level, cell.neighbors);
  fclose (fp);

  sprintf (name, "ghost-%d", pid());
  fp = fopen (name, "w");
  for (int i = 0; i < m->nsnd; i++)
    for (int j = 0; j <= m->snd[i].depth; j++)
      foreach_cache_level (m->snd[i].halo[j], j,)
	fprintf (fp, "%g %g %d %d\n", x, y, m->snd[i].pid,
		 allocated(0,0) ? cell.neighbors : -1);
  fclose (fp);
#endif
}

void mpi_boundary_update (void * p)
{
  MpiBoundary * m = p ? p : (MpiBoundary *) mpi_boundary;

  /* send refinement status of halo cells */
  int * buf[m->nsnd];
  MPI_Request r[m->nsnd];
  for (int i = 0; i < m->nsnd; i++) {
    Rcv snd = m->snd[i];
    size_t n = 0;
    for (int l = snd.depth; l >= 0; l--)
      n += snd.halo[l].n;
    assert (n > 0);
    buf[i] = malloc (sizeof (int)*n);
    int * b = buf[i];
    for (int l = snd.depth; l >= 0; l--)
      foreach_cache_level(snd.halo[l], l,)
	*b++ = !is_leaf(cell);
    MPI_Isend (buf[i], n, MPI_INT, snd.pid, 0, MPI_COMM_WORLD, &r[i]);
  }
    
  /* receive refinement status of halo cells */
  for (int i = 0; i < m->nrcv; i++) {
    Rcv * rcv = &m->rcv[i];
    size_t n = 0;
    for (int l = rcv->depth; l >= 0; l--)
      n += rcv->halo[l].n;
    assert (n > 0);
    int buf[n], * b = buf;
    MPI_Status s;
    MPI_Recv (buf, n, MPI_INT, rcv->pid, 0, MPI_COMM_WORLD, &s);
    for (int l = rcv->depth; l >= 0; l--)
      foreach_cache_level(rcv->halo[l], l,) {
	bool refined = *b++;
	// add newly-refined halo cells to halo buffer
	if (is_leaf(cell) && refined) {
	  cell.flags |= halo;
	  cell.flags &= ~leaf;
	  if (cell.neighbors) {
	    foreach_child() {
	      rcv_append (point, rcv);
	      cell.flags &= ~halo;
	      cell.flags |= leaf;
	      assert (!is_active(cell));
	    }
	    ((Quadtree *)grid)->dirty = true;
	  }
	}
      }
  }
    
  /* check that ghost values were received OK and free send buffers */
  for (int i = 0; i < m->nsnd; i++) {
    MPI_Status s;
    MPI_Wait (&r[i], &s);
    free (buf[i]);
  }
  
  mpi_boundary_sync (m);
}

void mpi_partitioning()
{
  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  int nf = 0;
  foreach()
    nf++;
  nf = max(1, nf/m->npe);

  // set the pid of each cell
  scalar pid[];
  int i = 0;
  foreach_cell_post (is_active (cell) || cell.neighbors > 0)
    if (is_active (cell)) {
      if (is_leaf (cell)) {
	pid[] = min(m->npe - 1, i/nf);
	i++;
      }
      else {
	pid[] = -1;
	for (int i = 0; i <= 1; i++)
	  for (int j = 0; j <= 1; j++)
	    if (pid[] == -1)
	      pid[] = fine(pid,i,j);
	    else if (fine(pid,i,j) == pid())
	      pid[] = pid();
      }
      if (pid[] != pid())
	cell.flags &= ~active;
    }

  /* this is the adjacency vector i.e. adjacency_rcv[x] is different from
     zero if we need to receive data from proc x */
  int adjacency_rcv[m->npe];
  for (int i = 0; i < m->npe; i++)
    adjacency_rcv[i] = 0;

  /* we build arrays of ghost cell indices for restriction */
  ((Quadtree *)grid)->dirty = true;
  foreach_cell() {
    if (!is_active(cell)) {
      if (!is_leaf(cell))
	// decrement_neighbors (point);
	for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	  for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	    if (allocated(k,l) && !is_active(neighbor(k,l))) {
	      neighbor(k,l).neighbors--;
	      if (neighbor(k,l).neighbors == 0)
		free_children (point, k, l);
	    }
      int nactive = 0;
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  if (allocated(o,p) && is_active(neighbor(o,p)))
	    nactive++;
      if (nactive > 0) {
	int i;
	for (i = 0; i < m->nrcv; i++)
	  if (pid[] == m->rcv[i].pid)
	    break;
	if (i == m->nrcv) {
	  m->rcv = realloc (m->rcv, ++m->nrcv*sizeof (Rcv));
	  Rcv * rcv = &m->rcv[m->nrcv-1];
	  rcv->pid = pid[];
	  rcv->depth = 0;
	  rcv->halo = malloc (sizeof (CacheLevel));
	  cache_level_init (&rcv->halo[0]);
	}
	rcv_append (point, &m->rcv[i]);
	adjacency_rcv[(int)pid[]]++;
      }
    }
    if (is_leaf(cell))
      continue;
  }

  /* we do a global transpose of the adjacency vector
     i.e. adjacency[x] is different from zero if we need to send
     data to proc x */
  MPI_Alltoall (adjacency_rcv, 1, MPI_INT, m->adjacency, 1, MPI_INT,
		MPI_COMM_WORLD);

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

