#define DEBUG 0

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
  void * buf;        // MPI buffer
  MPI_Request r;     // MPI request
} Rcv;

typedef struct {
  Rcv * rcv;
  int npid;
  int pid[sq(2*GHOSTS+1)], n;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;
  
  SndRcv restriction, halo_restriction, prolongation;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

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

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      foreach_cache_level(rcv->halo[i], i,)
	fprintf (fp, "%s%g %g %d %d\n", prefix, x, y, rcv->pid, level);
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("snd_rcv_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    free (rcv->buf);
    rcv->buf = NULL;
    prof_stop();
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

static RcvPid * rcv_pid_new (void)
{
  return calloc (sizeof (RcvPid), 1);
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());
  for (int i = 0; i < p->n; i++)
    if (p->pid[i] == pid)
      return NULL;

  assert (p->n < sq(2*GHOSTS+1));
  p->pid[p->n++] = pid;
  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = realloc (p->rcv, ++p->npid*sizeof (Rcv));
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = 0;
    rcv->halo = malloc (sizeof (CacheLevel));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, pid);
  if (rcv)
    rcv_append (point, rcv);
}

static void rcv_pid_append_children (RcvPid * p, int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, pid);
  if (rcv)
    foreach_child()
      rcv_append (point, rcv);
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  free (p->rcv);
  free (p);
}

static Boundary * mpi_boundary = NULL;

static void rcv_pid_receive (RcvPid * m, scalar * list, int l)
{
  int len = list_len (list);
    
  /* receive ghost values */
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      prof_start ("rcv_pid_receive");
      double buf[rcv->halo[l].n*len];
#if 0
      fprintf (stderr, "receiving %d doubles from %d level %d\n",
	       rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      MPI_Recv (buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid, l, 
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      double * b = buf;
      foreach_cache_level(rcv->halo[l], l,)
	for (scalar s in list)
	  s[] = *b++;
      prof_stop();
    }
  }
}

static void rcv_pid_wait (RcvPid * m, int l)
{
  /* wait for completion of send requests */
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      rcv_free_buf (rcv);
  }
}

static void rcv_pid_send (RcvPid * m, scalar * list, int l)
{
  int len = list_len (list); 

  /* send ghost values */
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      prof_start ("rcv_pid_send");
      assert (!rcv->buf);
      rcv->buf = malloc (sizeof (double)*rcv->halo[l].n*len);
      double * b = rcv->buf;
      foreach_cache_level(rcv->halo[l], l,)
	for (scalar s in list)
	  *b++ = s[];
#if 0
      fprintf (stderr, "sending %d doubles to %d level %d\n",
	       rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      MPI_Isend (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid, l, 
		 MPI_COMM_WORLD, &rcv->r);
      prof_stop();
    }
  }
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  rcv_pid_send (m->snd, list, l);
  rcv_pid_receive (m->rcv, list, l);
  rcv_pid_wait (m->snd, l);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m)
{
  m->rcv = rcv_pid_new();
  m->snd = rcv_pid_new();
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->restriction);
  snd_rcv_destroy (&m->halo_restriction);
  snd_rcv_destroy (&m->prolongation);
  free (m);
}

static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
}

static void mpi_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->halo_restriction, list, l);
}

static void mpi_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, int l, int depth)
{
  MpiBoundary * m = (MpiBoundary *) b;
  if (l == depth)
    rcv_pid_sync (&m->restriction, list, l);
  else
    rcv_pid_sync (&m->prolongation, list, l);
}

static void mpi_boundary_halo_restriction_flux (const Boundary * b,
						vector * list)
{
  //  MpiBoundary * mpi = (MpiBoundary *) b;
}

void mpi_boundary_new()
{
  mpi_boundary = calloc (1, sizeof (MpiBoundary));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->restriction = mpi_boundary_restriction;
  mpi_boundary->halo_restriction = mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation = mpi_boundary_halo_prolongation;
  mpi_boundary->halo_restriction_flux = mpi_boundary_halo_restriction_flux;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->restriction);
  snd_rcv_init (&mpi->halo_restriction);
  snd_rcv_init (&mpi->prolongation);
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells (FILE * fp);

  char prefix[80];
  FILE * fp;
  
  // local halo
  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l)
      foreach_child()
        fprintf (fp, "%s%g %g %d\n", prefix, x, y, level);
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }
  
  fp = fopen_prefix (fp1, "neighbors", prefix);
  foreach()
    fprintf (fp, "%s%g %g %d\n", prefix, x, y, cell.neighbors);
  if (!fp1)
    fclose (fp);

  // local restriction
  fp = fopen_prefix (fp1, "restriction", prefix);
  for (int l = 0; l < depth(); l++)
    foreach_halo (restriction, l)
      fprintf (fp, "%s%g %g %d %d\n", prefix, x, y, level, cell.neighbors);
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  
  fp = fopen_prefix (fp1, "mpi-restriction", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-halo-restriction", prefix);
  rcv_pid_print (mpi->halo_restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-prolongation", prefix);
  rcv_pid_print (mpi->prolongation.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-border", prefix);
  foreach_cell() {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %d %d %d\n",
	       prefix, x, y, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
  foreach_cell() {
    if (!is_active(cell))
      fprintf (fp, "%s%g %g %d %d %d\n",
	       prefix, x, y, level, cell.neighbors, cell.pid);
    else if (!is_border(cell))
      continue;
    if (is_leaf(cell))
      continue;
  }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new();
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new();
}

void mpi_boundary_update (MpiBoundary * m)
{
  prof_start ("mpi_boundary_update");

  SndRcv * restriction = &m->restriction;
  SndRcv * prolongation = &m->prolongation;
  SndRcv * halo_restriction = &m->halo_restriction;

  snd_rcv_free (restriction);
  snd_rcv_free (prolongation);
  snd_rcv_free (halo_restriction);
  
  #define is_prolongation(cell) (!is_leaf(cell) && !cell.neighbors)

  /* we build arrays of ghost cell indices for restriction and prolongation */
  foreach_cell() {
    if (is_active(cell) && !is_border(cell))
      // we skip the interior of the local domain
      continue;

    if (is_local(cell)) {
      // ==================================
      // local cell: do we need to send it?
      RcvPid * pro = prolongation->snd;
      RcvPid * res = restriction->snd;
      pro->n = res->n = 0;
      for (int k = -GHOSTS; k <= GHOSTS; k++)
	for (int l = -GHOSTS; l <= GHOSTS; l++) {
	  int pid = neighbor(k,l).pid;
	  if (pid >= 0 && pid != cell.pid && !is_prolongation(neighbor(k,l))) {
	    rcv_pid_append (res, pid, point);
	    if (is_leaf(neighbor(k,l)))
	      rcv_pid_append (pro, pid, point);
	  }
	}
      // prolongation
      if (is_leaf(cell)) {
	if (cell.neighbors) {
	  pro->n = res->n = 0;
	  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	      if (neighbor(k,l).pid >= 0 &&
		  !is_leaf(neighbor(k,l)) && neighbor(k,l).neighbors)
		for (int i = 2*k; i <= 2*k + 1; i++)
		  for (int j = 2*l; j <= 2*l + 1; j++)
		    if (child(i,j).pid != cell.pid) {
		      rcv_pid_append_children (res, child(i,j).pid, point);
		      if (is_leaf(child(i,j)))
			rcv_pid_append_children (pro, child(i,j).pid, point);
		    }
	}
      }
      else
	// is it the parent of a non-local cell?
	for (int k = -2; k <= 3; k++)
	  for (int l = -2; l <= 3; l++) {
	    int pid = child(k,l).pid;
	    if (pid >= 0 && pid != cell.pid)
	      rcv_pid_append (res, pid, point);
	  }
      // halo restriction
      if (level > 0 && !is_local(aparent(0,0))) {
	RcvPid * halo_res = halo_restriction->snd;
	halo_res->n = 0;
	rcv_pid_append (halo_res, aparent(0,0).pid, point);
      }
    }
    else {
      // =========================================
      // non-local cell: do we need to receive it?
      RcvPid * pro = prolongation->rcv;
      RcvPid * res = restriction->rcv;
      pro->n = res->n = 0;
      for (int k = -GHOSTS; k <= GHOSTS; k++)
	for (int l = -GHOSTS; l <= GHOSTS; l++)
	  if (allocated(k,l) && is_active(neighbor(k,l)) &&
	      is_local(neighbor(k,l))) {
	    rcv_pid_append (res, cell.pid, point);
	    if (is_leaf(neighbor(k,l)))
	      rcv_pid_append (pro, cell.pid, point);
	  }
      if (is_leaf(cell)) {
	if (cell.neighbors) {
	  // prolongation
	  pro->n = res->n = 0;
	  for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	    for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	      if (allocated(k,l) && is_active(neighbor(k,l)) &&
		  !is_leaf(neighbor(k,l))) {
		for (int i = 2*k; i <= 2*k + 1; i++)
		  for (int j = 2*l; j <= 2*l + 1; j++)
		    if (is_local(child(i,j))) {
		      rcv_pid_append_children (res, cell.pid, point);
		      if (is_leaf(child(i,j)))
			rcv_pid_append_children (pro, cell.pid, point);
		    }
	      }
	}
      }
      else
	// is it the parent of a local cell?
	for (int k = -2; k <= 3; k++)
	  for (int l = -2; l <= 3; l++)
	    if (is_local(child(k,l)))
	      rcv_pid_append (res, cell.pid, point);
      // halo restriction
      if (level > 0 && is_local(aparent(0,0))) {
	RcvPid * halo_res = halo_restriction->rcv;
	halo_res->n = 0;
	rcv_pid_append (halo_res, cell.pid, point);
      }
    }

    if (is_leaf(cell))
      continue;
  }

  // fixme: cleanup cells which do not have local neighbors
  
  prof_stop();

#if DEBUG
  debug_mpi (NULL);
#endif  
}

/*
  Returns a linear quadtree containing the leaves which are in the
  (5x5) neighborhood of cells belonging to the (remote) process 'rpid'.
 */
Array * remote_leaves (unsigned rpid)
{
  static const int remote = 1 << user;
  foreach_cell_post (is_border(cell) && !is_leaf(cell))
    if (is_border(cell)) {
      if (!(cell.flags & remote))
	for (int k = -GHOSTS; k <= GHOSTS; k++)
	  for (int l = -GHOSTS; l <= GHOSTS; l++)
	    if (neighbor(k,l).pid == rpid) {
	      cell.flags |= remote;
	      k = l = GHOSTS + 1; // break
	    }
      if (level > 0 && (cell.flags & remote))
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
static int match_refine (unsigned * a, scalar * list)
{
  int nactive = 0;
  if (a != NULL)
    foreach_cell() {
      unsigned flags = *a++;
      if (flags & leaf)
	continue;
      else if (is_leaf(cell))
	refine_cell (point, list, 0, &nactive);	
    }
  return nactive;
}

static void mpi_boundary_match (MpiBoundary * mpi)
{
  prof_start ("mpi_boundary_match");
  
  /* Send halo mesh for each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  Array * a[snd->npid];
  MPI_Request r[snd->npid];
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid;
    a[i] = remote_leaves (pid);
    assert (a[i]->len > 0);
    MPI_Isend (a[i]->p, a[i]->len, MPI_INT, pid, 0, MPI_COMM_WORLD, &r[i]);
  }
    
  /* Receive halo mesh from each neighboring process. */
  
  RcvPid * rcv = mpi->restriction.rcv;
  int refined = 0;
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid;
    MPI_Status s;
    int len;
    MPI_Probe (pid, 0, MPI_COMM_WORLD, &s);
    MPI_Get_count (&s, MPI_INT, &len);
    unsigned a[len];
    MPI_Recv (a, len, MPI_INT, pid, 0, MPI_COMM_WORLD, &s);
    refined += match_refine (a, NULL);
  }
  
  /* check that ghost values were received OK and free send buffers */
  for (int i = 0; i < snd->npid; i++) {
    MPI_Wait (&r[i], MPI_STATUS_IGNORE);
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
}

void mpi_partitioning()
{
  prof_start ("mpi_partitioning");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  int nf = 0;
  foreach()
    nf++;
  nf = max(1, nf/npe());

  /* set the pid of each cell */
  int i = 0;
  quadtree->dirty = true;
  foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
	cell.pid = min(npe() - 1, i/nf); i++;
	if (cell.neighbors > 0)
	  for (int i = 0; i <= 1; i++)
	    for (int j = 0; j <= 1; j++)
	      child(i,j).pid = cell.pid;
	if (!is_local(cell))
	  cell.flags &= ~active;
      }
      else {
	cell.pid = child(0,1).pid;
	bool inactive = true;
	for (int i = 0; i <= 1; i++)
	  for (int j = 0; j <= 1; j++)
	    if (is_active(child(i,j)))
	      inactive = false;
	if (inactive)
	  cell.flags &= ~active;
      }
    }

  // flag border cells
  foreach_cell() {
    if (is_active(cell)) {
      for (int k = -GHOSTS; k <= GHOSTS; k++)
	for (int l = -GHOSTS; l <= GHOSTS; l++)
	  if (allocated(k,l) && !is_local(neighbor(k,l))) {
	    cell.flags |= border;
	    k = l = GHOSTS + 1;
	  }
    }
    else
      continue;
    if (is_leaf(cell)) {      
      if (is_border(cell) && cell.neighbors)
	for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	  for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	    if (allocated(k,l) && !is_local(neighbor(k,l))) {
	      for (int k = 0; k < 2; k++)
		for (int l = 0; l < 2; l++)
		  child(k,l).flags |= border;
	      k = l = GHOSTS + 1;
	    }
      continue;
    }
  }
  
  prof_stop();
  
  mpi_boundary_update (m);
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
    foreach_coarse_level(l) // fixme: we could use restriction()
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

