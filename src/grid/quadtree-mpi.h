#define DEBUG 0

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
  void * buf;        // MPI buffer
  MPI_Request r;     // MPI request
} Rcv;

#if dimension == 2
# define PIDMAX (sq(2*GHOSTS + 1) + sq(6))
#else
# define PIDMAX (cube(2*GHOSTS + 1) + cube(6))
#endif

typedef struct {
  int pid[PIDMAX], n;  
} PidArray;

typedef struct {
  Rcv * rcv;
  int npid;
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
	fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level);
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
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

static Rcv * rcv_pid_pointer (RcvPid * p, PidArray * a, int pid)
{
  assert (pid >= 0 && pid < npe());
  for (int i = 0; i < a->n; i++)
    if (a->pid[i] == pid)
      return NULL;

  assert (a->n < PIDMAX);
  a->pid[a->n++] = pid;
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

static void rcv_pid_append (RcvPid * p, PidArray * a,
			    int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, a, pid);
  if (rcv)
    rcv_append (point, rcv);
}

static void rcv_pid_append_children (RcvPid * p, PidArray * a,
				     int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, a, pid);
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

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*dimension*vectors_len (listf);

  /* initiate non-blocking receives */
  MPI_Request r[m->npid];
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    r[i] = MPI_REQUEST_NULL;
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = malloc (sizeof (double)*rcv->halo[l].n*len);
#if 0
      fprintf (stderr, "receiving %d doubles from %d level %d\n",
	       rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid, l,
		 MPI_COMM_WORLD, &r[i]);
    }
  }

  /* receive ghost values */
  int i;
  MPI_Status s;
  MPI_Waitany (m->npid, r, &i, &s);
  while (i != MPI_UNDEFINED) {
    Rcv * rcv = &m->rcv[i];
    assert (l <= rcv->depth && rcv->halo[l].n > 0);
    assert (rcv->buf);
    double * b = rcv->buf;
    foreach_cache_level(rcv->halo[l], l,) {
      for (scalar s in list)
	s[] = *b++;
      for (vector v in listf)
	foreach_dimension() {
	  if (allocated(-1) && !is_local(neighbor(-1)))
	    v.x[] = *b;
	  b++;
	  if (allocated(1) && !is_local(neighbor(1)))
	    v.x[1] = *b;
	  b++;
	}
    }
    int rlen;
    MPI_Get_count (&s, MPI_DOUBLE, &rlen);
    assert (rlen == (b - (double *) rcv->buf));
    free (rcv->buf);
    rcv->buf = NULL;
    MPI_Waitany (m->npid, r, &i, &s);
  }

  prof_stop();
}

static void rcv_pid_wait (RcvPid * m)
{
  /* wait for completion of send requests */
  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
}

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l)
{
  prof_start ("rcv_pid_send");

  int len = list_len (list) + 2*dimension*vectors_len (listf);

  /* send ghost values */
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = malloc (sizeof (double)*rcv->halo[l].n*len);
      double * b = rcv->buf;
      foreach_cache_level(rcv->halo[l], l,) {
	for (scalar s in list)
	  *b++ = s[];
	for (vector v in listf)
	  foreach_dimension() {
	    *b++ = v.x[];
	    *b++ = allocated(1) ? v.x[1] : undefined;
	  }
      }
#if 0
      fprintf (stderr, "sending %d doubles to %d level %d\n",
	       rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
		 MPI_DOUBLE, rcv->pid, l, MPI_COMM_WORLD, &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL;
  vector * listf = NULL;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face)
	listf = vectors_add (listf, s.v);
      else
	listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listf, l);
  rcv_pid_receive (m->rcv, listr, listf, l);
  rcv_pid_wait (m->snd);
  free (listr);
  free (listf);
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

trace
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
}

trace
static void mpi_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->halo_restriction, list, l);
}

trace
static void mpi_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, int l, int depth)
{
  MpiBoundary * m = (MpiBoundary *) b;
  if (l == depth)
    rcv_pid_sync (&m->restriction, list, l);
  else
    rcv_pid_sync (&m->prolongation, list, l);
}

void mpi_boundary_new()
{
  mpi_boundary = calloc (1, sizeof (MpiBoundary));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->restriction = mpi_boundary_restriction;
  mpi_boundary->halo_restriction = mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation = mpi_boundary_halo_prolongation;
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
      fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }
  
  fp = fopen_prefix (fp1, "neighbors", prefix);
  foreach()
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
  if (!fp1)
    fclose (fp);

  // local restriction
  fp = fopen_prefix (fp1, "restriction", prefix);
  for (int l = 0; l < depth(); l++)
    foreach_halo (restriction, l)
      fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, level, cell.neighbors);
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  
  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-halo-restriction-rcv", prefix);
  rcv_pid_print (mpi->halo_restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-prolongation-rcv", prefix);
  rcv_pid_print (mpi->prolongation.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);
  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-halo-restriction-snd", prefix);
  rcv_pid_print (mpi->halo_restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-prolongation-snd", prefix);
  rcv_pid_print (mpi->prolongation.snd, fp, prefix);
  if (!fp1)
    fclose (fp);
  
  fp = fopen_prefix (fp1, "mpi-border", prefix);
  foreach_cell() {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
	       prefix, x, y, z, level, cell.neighbors, cell.pid);
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
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
	       prefix, x, y, z, level, cell.neighbors, cell.pid,
	       is_remote_leaf(cell));
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

trace
void mpi_boundary_update()
{
  if (npe() == 1)
    return;

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  
  prof_start ("mpi_boundary_update");

  SndRcv * restriction = &m->restriction;
  SndRcv * prolongation = &m->prolongation;
  SndRcv * halo_restriction = &m->halo_restriction;

  snd_rcv_free (restriction);
  snd_rcv_free (prolongation);
  snd_rcv_free (halo_restriction);
  
  /* we build arrays of ghost cell indices for restriction and prolongation
     we do a breadth-first traversal from fine to coarse, so that
     coarsening of unused cells can proceed fully. See also
     quadtree-common.h:coarsen_function() */
  for (int l = depth(); l >= 0; l--)
    foreach_cell() {
      if (is_active(cell) && !is_border(cell))
	// we skip the interior of the local domain
	continue;

      if (level == l) {
	if (is_local(cell)) {
	  // ==================================
	  // local cell: do we need to send it?
	  RcvPid * pro = prolongation->snd;
	  RcvPid * res = restriction->snd;
	  PidArray apro, ares;
	  apro.n = ares.n = 0;
	  for (int k = -GHOSTS; k <= GHOSTS; k++)
	    for (int l = -GHOSTS; l <= GHOSTS; l++)
	    #if dimension == 3
	      for (int n = -GHOSTS; n <= GHOSTS; n++)
	    #endif
		{
		  int pid = neighbor(k,l,n).pid;
		  if (pid >= 0 && pid != cell.pid) {
		    rcv_pid_append (res, &ares, pid, point);
		    if ((is_leaf(neighbor(k,l,n)) ||
			 is_prolongation(neighbor(k,l,n))) ||
			(is_leaf(cell) &&
			 k >= -GHOSTS/2 && k <= GHOSTS/2 &&
			 l >= -GHOSTS/2 && l <= GHOSTS/2
			 #if dimension == 3
			 && n >= -GHOSTS/2 && n <= GHOSTS/2
			 #endif
			 ))
		      rcv_pid_append (pro, &apro, pid, point);
		  }
	      }
	  // prolongation
	  if (is_leaf(cell)) {
	    if (cell.neighbors) {
	      PidArray cpro, cres;
	      cpro.n = cres.n = 0;
	      for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
		for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
                #if dimension == 3
		  for (int n = -GHOSTS/2; n <= GHOSTS/2; n++)
                #endif
		    if (neighbor(k,l,n).pid >= 0 && neighbor(k,l,n).neighbors)
		      for (int i = 2*k; i <= 2*k + 1; i++)
			for (int j = 2*l; j <= 2*l + 1; j++)
			#if dimension == 3
			  for (int k = 2*n; k <= 2*n + 1; k++)
                        #endif
			    if (child(i,j,k).pid != cell.pid) {
			      rcv_pid_append_children (res, &cres,
						       child(i,j,k).pid, point);
			      if (!is_refined(child(i,j,k))) {
				rcv_pid_append (pro, &apro, child(i,j,k).pid,
						point);
				rcv_pid_append_children (pro, &cpro,
							 child(i,j,k).pid,
							 point);
			      }
			    }
	    }
	  }
	  else
	    // is it the parent of a non-local cell?
	    for (int k = -2; k <= 3; k++)
	      for (int l = -2; l <= 3; l++)
	      #if dimension == 3
		for (int n = -2; n <= 3; n++)
              #endif
		  {
		    int pid = child(k,l,n).pid;
		    if (pid >= 0 && pid != cell.pid)
		      rcv_pid_append (res, &ares, pid, point);
		  }
	  // halo restriction
	  if (level > 0 && !is_local(aparent(0))) {
	    RcvPid * halo_res = halo_restriction->snd;
	    apro.n = 0;
	    rcv_pid_append (halo_res, &apro, aparent(0).pid, point);
	  }
	}
	else {
	  // =========================================
	  // non-local cell: do we need to receive it?
	  RcvPid * pro = prolongation->rcv;
	  RcvPid * res = restriction->rcv;
	  bool neighbors = false;
	  PidArray apro, ares;
	  apro.n = ares.n = 0;
	  for (int k = -GHOSTS; k <= GHOSTS; k++)
	    for (int l = -GHOSTS; l <= GHOSTS; l++)
	    #if dimension == 3
	      for (int n = -GHOSTS; n <= GHOSTS; n++)
            #endif
		{
		  if (allocated(k,l,n) && is_local(neighbor(k,l,n))) {
		    neighbors = true;
		    rcv_pid_append (res, &ares, cell.pid, point);
		    if ((is_leaf(neighbor(k,l,n)) ||
			 !is_active(neighbor(k,l,n))) ||
			(is_leaf(cell) &&
			 k >= -GHOSTS/2 && k <= GHOSTS/2 &&
			 l >= -GHOSTS/2 && l <= GHOSTS/2
			 #if dimension == 3
			 && n >= -GHOSTS/2 && n <= GHOSTS/2
			 #endif
			 ))
		      rcv_pid_append (pro, &apro, cell.pid, point);
		  }
		}
	  if (is_leaf(cell)) {
	    if (cell.neighbors) {
	      // prolongation
	      PidArray cpro, cres;
	      cpro.n = cres.n = 0;
	      for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
		for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
		#if dimension == 3
		  for (int n = -GHOSTS/2; n <= GHOSTS/2; n++)
		#endif
		    if (allocated(k,l,n) && is_active(neighbor(k,l,n)) &&
			neighbor(k,l,n).neighbors)
		      for (int i = 2*k; i <= 2*k + 1; i++)
			for (int j = 2*l; j <= 2*l + 1; j++)
			#if dimension == 3
			  for (int k = 2*n; k <= 2*n + 1; k++)
                        #endif
			    if (is_local(child(i,j,k))) {
			      rcv_pid_append_children (res, &cres, cell.pid,
						       point);
			      if (!is_refined(child(i,j,k))) {
				rcv_pid_append (pro, &apro, cell.pid, point);
				rcv_pid_append_children (pro, &cpro,
							 cell.pid, point);
			      }
			    }
	    }
	  }
	  else {
	    // not a leaf
	    // is it the parent of a local cell?
	    for (int k = -2; k <= 3; k++)
	      for (int l = -2; l <= 3; l++)
	      #if dimension == 3
		for (int n = -2; n <= 3; n++)
              #endif
		  if (is_local(child(k,l,n))) {
		    neighbors = true;
		    rcv_pid_append (res, &ares, cell.pid, point);
		  }
	    // coarse cell with no neighbors: destroy its children
	    if (!neighbors)
	      coarsen_cell (point, NULL, NULL);
	  }
	  // halo restriction
	  if (level > 0 && is_local(aparent(0))) {
	    RcvPid * halo_res = halo_restriction->rcv;
	    apro.n = 0;
	    rcv_pid_append (halo_res, &apro, cell.pid, point);
	  }
	}
	continue; // level == l
      }
    
      if (is_leaf(cell))
	continue;
    }
  
  prof_stop();

#if DEBUG
  debug_mpi (NULL);
#endif  
}

/*
  Returns a linear quadtree containing the cells which are in the
  (5x5) neighborhood of cells belonging to the (remote) process 'rpid'.
 */
trace
Array * remote_cells (unsigned rpid)
{
  static const int remote = 1 << user;
  foreach_cell_post (is_border(cell) && !is_leaf(cell))
    if (is_border(cell)) {
      if (!(cell.flags & remote)) {
	short flags = cell.flags;
	foreach_neighbor()
	  if (cell.pid == rpid) {
	    flags |= remote;
	    break;
	  }
	cell.flags = flags;
      }
      if (level > 0 && (cell.flags & remote))
	aparent(0).flags |= remote;
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
    if (is_leaf(cell))
      flags |= remote_leaf;
    array_append (a, &flags);
    if (flags & leaf)
      continue;
  }
  return a;
}

/**
   Match the refinement given by the linear quadtree 'a'. 
   Returns the number of border cells refined.
*/
trace
int match_refine (unsigned * a, scalar * list, int pid)
{
  int nr = 0;
  if (a != NULL)
    foreach_cell() {
      unsigned flags = *a++;
      if (cell.pid == pid) {
	if (flags & remote_leaf)
	  cell.flags |= remote_leaf;
	else
	  cell.flags &= ~remote_leaf;
      }
      if (flags & leaf)
	continue;
      else if (is_leaf(cell))
	nr += refine_cell (point, list, 0, NULL);
    }
  return nr;
}

#if 0
trace
void mpi_boundary_refine (scalar * list)
{
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  
  /* Send halo mesh for each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  Array * a[snd->npid];
  MPI_Request r[2*snd->npid];
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid;
    a[i] = remote_cells (pid);
    assert (a[i]->len > 0);
    int len = a[i]->len;
    MPI_Isend (&len, 1, MPI_INT, pid, 0,  MPI_COMM_WORLD, &r[2*i]);
    MPI_Isend (a[i]->p, len, MPI_INT, pid, 0, MPI_COMM_WORLD, &r[2*i+1]);
  }
    
  /* Receive halo mesh from each neighboring process. */
  RcvPid * rcv = mpi->restriction.rcv;
  int refined = 0;
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid;
    // we could use MPI_Probe, but apparently it's not a good idea,
    // see http://cw.squyres.com/ and
    // http://cw.squyres.com/columns/2004-07-CW-MPI-Mechanic.pdf
    int len;
    MPI_Recv (&len, 1, MPI_INT, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    unsigned a[len];
    MPI_Recv (a, len, MPI_INT, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    refined += match_refine (a, list, pid);
  }
  
  /* check that ghost values were received OK and free send buffers */
  MPI_Waitall (2*snd->npid, r, MPI_STATUSES_IGNORE);
  for (int i = 0; i < snd->npid; i++)
    array_free (a[i]);
  
  prof_stop();
  
  /* If an active cell has been refined (due to the 2:1 refinement
     constraint), we need to repeat the process. */
  mpi_all_reduce (refined, MPI_INT, MPI_SUM);
  if (refined)
    mpi_boundary_refine (list);
}
#else
trace
void mpi_boundary_refine (scalar * list)
{
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  /* Send refinement cache to each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  MPI_Request r[2*snd->npid];
  int nr = 0;
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid;
    int len = quadtree->refined.n;
    MPI_Isend (&len, 1, MPI_INT, pid, 0,  MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (quadtree->refined.p, sizeof(Index)/sizeof(int)*len,
		 MPI_INT, pid, 0, MPI_COMM_WORLD, &r[nr++]);
  }

  /* Receive refinement cache from each neighboring process. 
   fixme: use non-blocking receives */
  RcvPid * rcv = mpi->restriction.rcv;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid, len;
    MPI_Recv (&len, 1, MPI_INT, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (len > 0) {
      Index p[len];
      MPI_Recv (p, sizeof(Index)/sizeof(int)*len,
		MPI_INT, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      Cache refined = {p, len, len};
      foreach_cache (refined,)
	if (allocated(0)) {
	  if (is_leaf(cell)) {
	    bool neighbors = false;
	    foreach_neighbor()
	      if (allocated(0) && is_local(cell)) {
		neighbors = true;
		break;
	      }
	    // is this a non-local prolongation?
	    if (!neighbors)
	      foreach_neighbor (GHOSTS/2) {
		if (allocated(0) && is_active(cell) && cell.neighbors)
		  foreach_child()
		    if (is_local(cell)) {
		      neighbors = true;
		      break;
		    }
		if (neighbors)
		  break;
	      }
	    // refine the cell only if it has local neighbors
	    if (neighbors)
	      refine_cell (point, list, 0, &rerefined);
	  }
	  if (is_remote_leaf(cell))
	    cell.flags &= ~remote_leaf;
	}
    }
  }
  
  /* check that caches were received OK and free ressources */
  MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);

  /* update the refinement cache with "re-refined" cells */
  free (quadtree->refined.p);
  quadtree->refined = rerefined;
  
  prof_stop();

  /* if any cell has been re-refined, we repeat the process to take care
     of recursive refinements induced by the 2:1 constraint */
  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
}
#endif

trace
void mpi_boundary_coarsen (int l)
{
  prof_start ("mpi_boundary_coarsen");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  
  /* Send coarsening cache to each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  MPI_Request r[2*snd->npid];
  int nr = 0;
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid;
    int len = quadtree->coarsened.n;
    MPI_Isend (&len, 1, MPI_INT, pid, l,  MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (quadtree->coarsened.p, 2*len, MPI_INT, pid, l,
		 MPI_COMM_WORLD, &r[nr++]);
  }

  /* Receive coarsening cache from each neighboring process. 
   fixme: use non-blocking receives */
  RcvPid * rcv = mpi->restriction.rcv;
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid, len;
    MPI_Recv (&len, 1, MPI_INT, pid, l, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (len > 0) {
      IndexLevel p[len];
      MPI_Recv (p, 2*len, MPI_INT, pid, l, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      CacheLevel coarsened = {p, len, len};
      foreach_cache_level (coarsened, l,)
	if (allocated(0)) {
	  if (is_refined(cell))
	    assert (coarsen_cell (point, NULL, NULL));
	  if (cell.pid == pid && is_leaf(cell) && !is_remote_leaf(cell))
	    cell.flags |= remote_leaf;
	}
    }
  }
  
  /* check that caches were received OK and free ressources */
  MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);
  
  prof_stop();  
}

trace
void mpi_partitioning()
{
  prof_start ("mpi_partitioning");

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
	if (cell.neighbors > 0) {
	  int pid = cell.pid;
	  foreach_child()
	    cell.pid = pid;
	}
	if (!is_local(cell)) {
	  cell.flags &= ~active;
	  cell.flags |= remote_leaf;
	}
      }
      else {
	cell.pid = child(0,1,1).pid;
	bool inactive = true;
	foreach_child()
	  if (is_active(cell)) {
	    inactive = false;
	    break;
	  }
	if (inactive)
	  cell.flags &= ~(active|remote_leaf);
      }
    }

  // flag border cells
  foreach_cell() {
    if (is_active(cell)) {
      bool notb = true;
      for (int k = -GHOSTS; k <= GHOSTS && notb; k++)
	for (int l = -GHOSTS; l <= GHOSTS && notb; l++)
	#if dimension == 3
	  for (int n = -GHOSTS; n <= GHOSTS && notb; n++)
	#endif
	    if (allocated(k,l,n) && !is_local(neighbor(k,l,n))) {
	      cell.flags |= border;
	      notb = false;
	    }
    }
    else
      continue;
    if (is_leaf(cell)) {      
      if (is_border(cell) && cell.neighbors) {
	bool notb = true;
	for (int k = -GHOSTS/2; k <= GHOSTS/2 && notb; k++)
	  for (int l = -GHOSTS/2; l <= GHOSTS/2 && notb; l++)
	  #if dimension == 3
	    for (int n = -GHOSTS/2; n <= GHOSTS/2 && notb; n++)
          #endif
	      if (allocated(k,l,n) && !is_local(neighbor(k,l,n))) {
		foreach_child()
		  cell.flags |= border;
		notb = false;
	      }
      }
      continue;
    }
  }
  
  prof_stop();
  
  mpi_boundary_update();
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

trace
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
    foreach_coarse_level(l) { // fixme: we could use restriction()
      double sum = 0.;
      foreach_child()
	sum += size[];
      size[] = sum;
    }
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

