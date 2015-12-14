#define DEBUG 0

// this can be used to control the outputs of debug_mpi()
int debug_iteration = -1;

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  void * buf;        // MPI buffer
  MPI_Request r;     // MPI request
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
  int maxdepth;      // the maximum depth for this PE (= depth or depth + 1)
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
  char * name;
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
  if (is_refined(cell))
    level++;
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
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

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = calloc (sizeof (RcvPid), 1);
  r->name = strdup (name);
  return r;
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
    rcv->depth = rcv->maxdepth = 0;
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
  free (p->name);
  free (p);
}

static Boundary * mpi_boundary = NULL;

#define BOUNDARY_TAG(level) (level)
#define COARSEN_TAG(level)  ((level) + 64)
#define REFINE_TAG()        (128)

static size_t apply_bc (Rcv * rcv, scalar * list, vector * listf, int l)
{
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
  size_t size = b - (double *) rcv->buf;
  free (rcv->buf);
  rcv->buf = NULL;
  return size;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*dimension*vectors_len (listf);

  MPI_Request r[m->npid];
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    r[i] = MPI_REQUEST_NULL;
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = malloc (sizeof (double)*rcv->halo[l].n*len);
#if 0
      fprintf (stderr, "%s receiving %d doubles from %d level %d\n",
	       m->name, rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
#if 1 /* initiate non-blocking receive */
      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
		 BOUNDARY_TAG(l), MPI_COMM_WORLD, &r[i]);
#else /* blocking receive (useful for debugging) */
      MPI_Recv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
		BOUNDARY_TAG(l), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      apply_bc (rcv, list, listf, l);
#endif
    }
  }

  /* non-blocking receives (does nothing when using blocking receives) */
  int i;
  MPI_Status s;
  MPI_Waitany (m->npid, r, &i, &s);
  while (i != MPI_UNDEFINED) {
    Rcv * rcv = &m->rcv[i];
    assert (l <= rcv->depth && rcv->halo[l].n > 0);
    assert (rcv->buf);
    size_t size = apply_bc (rcv, list, listf, l);
    int rlen;
    MPI_Get_count (&s, MPI_DOUBLE, &rlen);
    assert (rlen == size);
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
      fprintf (stderr, "%s sending %d doubles to %d level %d\n",
	       m->name, rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
		 MPI_DOUBLE, rcv->pid, BOUNDARY_TAG(l), MPI_COMM_WORLD, &rcv->r);
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

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
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
  snd_rcv_init (&mpi->restriction, "restriction");
  snd_rcv_init (&mpi->halo_restriction, "halo_restriction");
  snd_rcv_init (&mpi->prolongation, "prolongation");
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
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
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

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= restriction.snd ======\n");
  RcvPid * snd = mpi->restriction.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= restriction.rcv ======\n");
  snd = mpi->restriction.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
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
	  int pid = cell.pid;
	  Point p = point;
	  foreach_neighbor()
	    if (cell.pid >= 0 && cell.pid != pid) {
	      rcv_pid_append (res, &ares, cell.pid, p);
	      if (is_leaf(cell) || is_prolongation(cell))
		rcv_pid_append (pro, &apro, cell.pid, p);
	    }
	  if (is_leaf(cell)) {
	    foreach_neighbor (GHOSTS/2)
	      if (cell.pid >= 0 && cell.pid != pid)
		rcv_pid_append (pro, &apro, cell.pid, p);
	    if (cell.neighbors) {
	      // prolongation
	      PidArray cpro, cres;
	      cpro.n = cres.n = 0;
	      int pid = cell.pid;
	      Point p = point;
	      foreach_neighbor (GHOSTS/2)
		if (cell.pid >= 0 && cell.neighbors)
		  foreach_child()
		    if (cell.pid != pid) {
		      rcv_pid_append_children (res, &cres, cell.pid, p);
		      if (!is_refined(cell)) {
			rcv_pid_append (pro, &apro, cell.pid, p);
			rcv_pid_append_children (pro, &cpro, cell.pid, p);
		      }
		    }
	    }
	  }
	  else
	    // not a leaf
	    // is it the parent of a non-local cell?
	    foreach_neighbor (GHOSTS/2)
	      foreach_child()
	        if (cell.pid >= 0 && cell.pid != pid)
		  rcv_pid_append (res, &ares, cell.pid, p);
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
	  int pid = cell.pid;
	  Point p = point;
	  foreach_neighbor()
	    if (allocated(0) && is_local(cell)) {
	      neighbors = true;
	      rcv_pid_append (res, &ares, pid, p);
	      if (is_leaf(cell) || !is_active(cell))
		rcv_pid_append (pro, &apro, pid, p);
	    }
	  if (is_leaf(cell)) {
	    foreach_neighbor (GHOSTS/2)
	      if (allocated(0) && is_local(cell)) {
		rcv_pid_append (pro, &apro, pid, p);
		break;
	      }
	    if (cell.neighbors) {
	      // prolongation
	      PidArray cpro, cres;
	      cpro.n = cres.n = 0;
	      int pid = cell.pid;
	      Point p = point;
	      foreach_neighbor (GHOSTS/2)
		if (allocated(0) && is_active(cell) && cell.neighbors)
		  foreach_child()
		    if (is_local(cell)) {
		      rcv_pid_append_children (res, &cres, pid, p);
		      if (!is_refined(cell)) {
			rcv_pid_append (pro, &apro, pid, p);
			rcv_pid_append_children (pro, &cpro, pid, p);
			break;
		      }
		    }
	    }
	  }
	  else {
	    // not a leaf
	    // is it the parent of a local cell?
	    foreach_neighbor (GHOSTS/2)
	      foreach_child()
	        if (is_local(cell)) {
		  neighbors = true;
		  rcv_pid_append (res, &ares, pid, p);
		  break;
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
    MPI_Isend (&len, 1, MPI_INT, pid, REFINE_TAG(),  MPI_COMM_WORLD, &r[2*i]);
    MPI_Isend (a[i]->p, len, MPI_INT, pid, REFINE_TAG(),
	       MPI_COMM_WORLD, &r[2*i+1]);
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
    MPI_Recv (&len, 1, MPI_INT, pid, REFINE_TAG(),
	      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    unsigned a[len];
    MPI_Recv (a, len, MPI_INT, pid, REFINE_TAG(),
	      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
    MPI_Isend (&len, 1, MPI_INT, pid, REFINE_TAG(),  MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (quadtree->refined.p, sizeof(Index)/sizeof(int)*len,
		 MPI_INT, pid, REFINE_TAG(), MPI_COMM_WORLD, &r[nr++]);
  }

  /* Receive refinement cache from each neighboring process. 
   fixme: use non-blocking receives */
  RcvPid * rcv = mpi->restriction.rcv;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0; i < rcv->npid; i++) {
    int pid = rcv->rcv[i].pid, len;
    MPI_Recv (&len, 1, MPI_INT, pid, REFINE_TAG(),
	      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (len > 0) {
      Index p[len];
      MPI_Recv (p, sizeof(Index)/sizeof(int)*len,
		MPI_INT, pid, REFINE_TAG(), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      Cache refined = {p, len, len};
      foreach_cache (refined,)
	if (level <= depth() && allocated(0)) {
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
  for (int i = 0; i < snd->npid; i++)
    if (snd->rcv[i].maxdepth > l) {
      int pid = snd->rcv[i].pid, len = quadtree->coarsened.n;
      MPI_Isend (&len, 1, MPI_INT, pid, COARSEN_TAG(l),
		 MPI_COMM_WORLD, &r[nr++]);
      if (len > 0)
	MPI_Isend (quadtree->coarsened.p, sizeof(IndexLevel)/sizeof(int)*len,
		   MPI_INT, pid, COARSEN_TAG(l), MPI_COMM_WORLD, &r[nr++]);
    }
    
  /* Receive coarsening cache from each (fine enough) neighboring process. 
   fixme: use non-blocking receives */
  RcvPid * rcv = mpi->restriction.rcv;
  for (int i = 0; i < rcv->npid; i++)
    if (rcv->rcv[i].maxdepth > l) {
      int pid = rcv->rcv[i].pid, len;
      MPI_Recv (&len, 1, MPI_INT, pid, COARSEN_TAG(l),
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (len > 0) {
	IndexLevel p[len];
	MPI_Recv (p, sizeof(IndexLevel)/sizeof(int)*len,
		  MPI_INT, pid, COARSEN_TAG(l),
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	CacheLevel coarsened = {p, len, len};
	foreach_cache_level (coarsened, l,) {
	  if (allocated(0)) {
	    if (is_refined(cell))
	      assert (coarsen_cell (point, NULL, NULL));
	    if (cell.pid == pid && is_leaf(cell) && !is_remote_leaf(cell))
	      cell.flags |= remote_leaf;
	  }
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
      short flags = cell.flags;
      foreach_neighbor()
	if (!is_local(cell)) {
	  flags |= border;
	  break;
	}
      cell.flags = flags;
    }
    else
      continue;
    if (is_leaf(cell)) {      
      if (is_border(cell) && cell.neighbors) {
	bool remote = false;
	foreach_neighbor (GHOSTS/2)
	  if (!is_local(cell)) {
	    remote = true;
	    break;
	  }
	if (remote)
	  foreach_child()
	    cell.flags |= border;
      }
      continue;
    }
  }
  
  prof_stop();
  
  mpi_boundary_update();
}

/**
# *z_indexing()*: fills *index* with the Z-ordering index.
   
If `leaves` is `true` only leaves are indexed, otherwise all active
cells are indexed. On a single processor, we would just need something
like (for leaves)

~~~literatec
double i = 0;
foreach()
  index[] = i++;
~~~

In parallel, this is a bit more difficult. 

On the master process (`pid() == 0`), the function returns the
(global) maximum index (and -1 on all other processes). */

trace
double z_indexing (scalar index, bool leaves)
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
      double sum = !leaves;
      foreach_child()
	sum += size[];
      size[] = sum;
    }
    boundary_iterate (restriction, {size}, l);
  }

  /**
  The maximum index value is the size of the entire tree (i.e. the
  value of `size` in the root cell on the master process) minus
  one. */
  
  double maxi = -1.;
  if (pid() == 0)
    foreach_level(0)
      maxi = size[] - 1.;
  
  /**
  ## Indexing
  
  Indexing can then be done locally. */

  double i = 0;
  foreach_cell() {
    double s = size[];
    if (!leaves || is_leaf(cell))
      index[] = i++;
    if (is_leaf(cell))
      continue;
    if (!is_active(cell)) {
      i += s - !leaves;
      continue;
    }
  }

  return maxi;
}
