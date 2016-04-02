// this can be used to control the outputs of debug_mpi()
int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo; // ghost cell indices for each level
  void * buf;        // MPI buffer
  MPI_Request r;     // MPI request
  int depth;         // the maximum number of levels
  int pid;           // the rank of the PE  
  int maxdepth;      // the maximum depth for this PE (= depth or depth + 1)
} Rcv;

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
  
  SndRcv restriction, restriction_root, halo_restriction, prolongation;
  Array * send, * receive; // which pids do we send to/receive from
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
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
      foreach_cache_level(rcv->halo[l], l,)
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

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());
  
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

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{
  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{
  // appends the pids of @p to @pids without duplication
  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
	break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
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
#define MOVED_TAG()         (256)

void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
		      MPI_Status s)
{
  double * b = rcv->buf;
  foreach_cache_level(rcv->halo[l], l,) {
    for (scalar s in list)
      s[] = *b++;
    for (vector v in listf) {
      foreach_dimension() {
	v.x[] = *b++;
	if (allocated(1))
	  v.x[1] = *b;
	b++;
      }
    }
  }
  size_t size = b - (double *) rcv->buf;
  free (rcv->buf);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (stderr,
	     "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
	     "Calling debug_mpi(NULL)...\n"
	     "Aborting...\n",
	     rlen, size, pid(), rcv->pid, l);
    fflush (stderr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}

#ifdef TIMEOUT
static struct {
  int count, source, tag;
  const char * name;
} RecvArgs;

static void timeout (int sig, siginfo_t *si, void *uc)
{
  fprintf (stderr,
	   "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n"
	   "Time out\n"
	   "Calling debug_mpi(NULL)...\n"
	   "Aborting...\n",
	   RecvArgs.name, RecvArgs.count, RecvArgs.source, RecvArgs.tag);
  fflush (stderr);
  debug_mpi (NULL);
  MPI_Abort (MPI_COMM_WORLD, -1);
}
#endif // TIMEOUT

static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
			    int source, int tag,
			    MPI_Comm comm, MPI_Status * status,
			    const char * name)
{
#ifdef TIMEOUT
  RecvArgs.count = count;
  RecvArgs.source = source;
  RecvArgs.tag = tag;
  RecvArgs.name = name;
  
  timer_t timerid;
  extern double t;
  if (t > TIMEOUT) {
    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = timeout;
    sigemptyset(&sa.sa_mask);
    assert (sigaction(SIGRTMIN, &sa, NULL) != -1);

    struct sigevent sev;
    sev.sigev_notify = SIGEV_SIGNAL;
    sev.sigev_signo = SIGRTMIN;
    sev.sigev_value.sival_ptr = &timerid;
    assert (timer_create(CLOCK_REALTIME, &sev, &timerid) != -1);

    struct itimerspec its;
    its.it_value.tv_sec = 10;
    its.it_value.tv_nsec = 0;
    its.it_interval.tv_sec = 0;
    its.it_interval.tv_nsec = 0;
    assert (timer_settime(timerid, 0, &its, NULL) != -1);
  }
#endif // TIMEOUT
  
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (stderr,
	     "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
	     "Calling debug_mpi(NULL)...\n"
	     "Aborting...\n",
	     name, count, source, tag, string);
    fflush (stderr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }

#ifdef TIMEOUT  
  if (t > TIMEOUT)
    assert (timer_delete (timerid) != -1);
#endif
}
  
static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;
  
  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*dimension*vectors_len (listf);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid]; // fixme: using NULL requests should be OK
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
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
		 BOUNDARY_TAG(l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;
#else /* blocking receive (useful for debugging) */
      MPI_Status s;
      mpi_recv_check (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
		      BOUNDARY_TAG(l), MPI_COMM_WORLD, &s, "rcv_pid_receive");
      apply_bc (rcv, list, listf, l, s);
#endif
    }
  }

  /* non-blocking receives (does nothing when using blocking receives) */
  if (nr > 0) {
    int i;
    MPI_Status s;
    MPI_Waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc (rcv, list, listf, l, s);
      MPI_Waitany (nr, r, &i, &s);
    }
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
  if (m->npid == 0)
    return;

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
		 MPI_DOUBLE, rcv->pid, BOUNDARY_TAG(l), MPI_COMM_WORLD,
		 &rcv->r);
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
  snd_rcv_destroy (&m->restriction_root);
  snd_rcv_destroy (&m->halo_restriction);
  snd_rcv_destroy (&m->prolongation);
  array_free (m->send);
  array_free (m->receive);
  free (m);
}

trace
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
  rcv_pid_sync (&m->restriction_root, list, l);
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
  if (true /*l == depth*/) {
    rcv_pid_sync (&m->restriction, list, l);
    rcv_pid_sync (&m->restriction_root, list, l);
  }
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
  snd_rcv_init (&mpi->restriction_root, "restriction_root");
  snd_rcv_init (&mpi->halo_restriction, "halo_restriction");
  snd_rcv_init (&mpi->prolongation, "prolongation");
  mpi->send = array_new();
  mpi->receive = array_new();
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

  // cleanup
  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "restriction-%d", pid()); remove (name);
    sprintf (name, "restriction-root-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-halo-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-prolongation-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-halo-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-prolongation-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }
  
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
  
  fp = fopen_prefix (fp1, "faces", prefix);
  foreach_face()
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
  foreach() {
    int n = 0;
    foreach_neighbor(1)
      if (is_refined(cell))
	n++;
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    assert (cell.neighbors == n);
  }
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

  fp = fopen_prefix (fp1, "mpi-restriction-root-rcv", prefix);
  rcv_pid_print (mpi->restriction_root.rcv, fp, prefix);
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

  fp = fopen_prefix (fp1, "mpi-restriction-root-snd", prefix);
  rcv_pid_print (mpi->restriction_root.snd, fp, prefix);
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
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
	       prefix, x, y, z, level, cell.neighbors,
	       cell.pid, cell.flags & leaf);
#if 0
    else if (is_active(cell) && !is_border(cell))
      continue;
    if (is_leaf(cell))
      continue;
#endif
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

  fp = fopen_prefix (fp1, "refined", prefix);
  foreach_cache (quadtree->refined,)
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);
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

static bool is_root (Point point)
{
  if (is_refined(cell))
    foreach_child()
      if (is_local(cell))
	return true;
  return false;
}

static void append_pid (Array * pids, Point point)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == cell.pid)
      return;
  array_append (pids, &cell.pid, sizeof(int));
}

// see src/figures/prolongation.svg
static bool is_local_prolongation (Point point, Point p)
{
#if dimension == 2
  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};
#elif dimension == 3
  struct { int x, y, z; } dp = {p.i - point.i, p.j - point.j, p.k - point.k};
#endif
  foreach_dimension() {
    if (dp.x == 0 && (is_refined(neighbor(-1)) || is_refined(neighbor(1))))
      return true;
    if (is_refined(neighbor(dp.x)))
      return true;
  }
  return false;
}

#define is_remote(cell) (cell.pid >= 0 && cell.pid != pid())

static int locals_pids (Point point, Array * pids)
{
  if (is_leaf(cell)) { // prolongation
    if (is_local(cell)) {
      Point p = point;
      foreach_neighbor(1) {
	if (is_remote(cell) &&
	    (is_refined(cell) || is_local_prolongation (point, p)))
	  append_pid (pids, point);
	if (is_refined(cell))
	  foreach_child()
	    if (is_remote(cell))
	      append_pid (pids, point);
      }
    }
  }
  else
    foreach_neighbor(1) {
      if (is_remote(cell))
	append_pid (pids, point);
      if (is_refined(cell))
	foreach_child()
	  if (is_remote(cell))
	    append_pid (pids, point);
    }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{
  foreach_child()
    if (is_remote(cell))
      append_pid (pids, point);
  return pids->len/sizeof(int);
}

// turns on quadtree_check() and co with outputs controlled by the condition
// #define DEBUGCOND (pid() >= 300 && pid() <= 400 && t > 0.0784876)

// turns on quadtree_check() and co without any outputs
// #define DEBUGCOND false

static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = malloc (npe()*sizeof(int));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) { // master
      // send/receive matrix i.e.
      // send[i][j] = number of points sent/received from i to j
      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
		  MPI_COMM_WORLD); 
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
		  MPI_COMM_WORLD);

      int * astatus = malloc (npe()*sizeof(int));
      for (int i = 0; i < npe(); i++)
	astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
	for (int j = 0; j < npe(); j++)
	  if (send[i][j] != receive[j][i]) {
	    fprintf (stderr, "%s: %d sends    %d to   %d at level %d\n",
		     name, i, send[i][j], j, l);
	    fprintf (stderr, "%s: %d receives %d from %d at level %d\n",
		     name, j, receive[j][i], i, l);
	    fflush (stderr);
	    for (int k = i - 2; k <= i + 2; k++)
	      if (k >= 0 && k < npe())
		astatus[k] = 1;
	    for (int k = j - 2; k <= j + 2; k++)
	      if (k >= 0 && k < npe())
		astatus[k] = 1;
	  }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      free (astatus);

      matrix_free (send);
      matrix_free (receive);
    }
    else { // slave
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD); 
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (stderr,
	       "check_snd_rcv_matrix \"%s\" failed\n"
	       "Calling debug_mpi(NULL)...\n"
	       "Aborting...\n",
	       name);
      fflush (stderr);
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  free (row);
}

trace
void mpi_boundary_update()
{
  if (npe() == 1)
    return;

  prof_start ("mpi_boundary_update");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * restriction = &m->restriction;
  SndRcv * restriction_root = &m->restriction_root;
  SndRcv * halo_restriction = &m->halo_restriction;

  snd_rcv_free (restriction);
  snd_rcv_free (restriction_root);
  snd_rcv_free (halo_restriction);
  
  static const unsigned short used = 1 << user;
  foreach_cell() {
    if (is_active(cell) && !is_border(cell))
      /* We skip the interior of the local domain.
	 Note that this can be commented out in case of suspicion that
	 something is wrong with border cell tagging. */
      continue;
    if (cell.neighbors) {
      // sending
      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
	foreach_child()
	  if (is_local(cell))
	    for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
	      rcv_pid_append (restriction->snd, *p, point);
	free (pids.p);
      }
      // receiving
      bool locals = false;
      if (is_leaf(cell)) { // prolongation
	if (is_remote(cell)) {
	  Point p = point;
	  foreach_neighbor(1)
	    if ((is_local(cell) &&
		 (is_refined(cell) || is_local_prolongation (point, p))) ||
		is_root(point))
	      locals = true, break;
	}
      }
      else
	foreach_neighbor(1)
	  if (is_local(cell) || is_root(point))
	    locals = true, break;
      if (locals)
	foreach_child()
	  if (is_remote(cell))
            rcv_pid_append (restriction->rcv, cell.pid, point),
	      cell.flags |= used;

      // root cells
      if (!is_leaf(cell)) {
	// sending
	if (is_local(cell)) {
	  Array pids = {NULL, 0, 0};
	  // root cell
	  int n = root_pids (point, &pids);
	  if (n) {
	    foreach_neighbor()
	      for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
		if (cell.pid >= 0 && cell.pid != *p)
		  rcv_pid_append (restriction_root->snd, *p, point);
	    free (pids.p);
	  }
	}
	// receiving
	else if (is_remote(cell)) {
	  bool root = false;
	  foreach_child()
	    if (is_local(cell))
	      root = true, break;
	  if (root) {
	    int pid = cell.pid;
	    foreach_neighbor()
	      if (is_remote(cell))
		rcv_pid_append (restriction_root->rcv, pid, point),
		  cell.flags |= used;
	  }
	}
      }
    }    
  }
    
  /* we remove unused cells
     we do a breadth-first traversal from fine to coarse, so that
     coarsening of unused cells can proceed fully. */
  
  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
    foreach_cell()
      if (level == l) {
	if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
	  aparent(0).flags |= keep;
	if (is_refined(cell) && !(cell.flags & keep))
	  coarsen_cell (point, NULL);
	cell.flags &= ~(used|keep);
	continue; // level == l
      }

  /* we update the list of send/receive pids */
  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (restriction->snd, m->send);
  rcv_pid_append_pids (restriction_root->snd, m->send);
  rcv_pid_append_pids (restriction->rcv, m->receive);
  rcv_pid_append_pids (restriction_root->rcv, m->receive);
  
  prof_stop();

#if DEBUG_MPI
  debug_mpi (NULL);
#endif

#ifdef DEBUGCOND
  extern double t;
  // check_snd_rcv_matrix (restriction, "restriction");
  // check_snd_rcv_matrix (restriction_root, "restriction_root");
  if (DEBUGCOND)
    debug_mpi (NULL);
  quadtree_check();
#endif
  
  // halo restriction
  scalar halov[];
  quadtree->dirty = true;
  for (int l = 0; l <= depth(); l++) {
    foreach_cell() {
      if (level == l) {
	if (is_local(cell)) {
	  bool restriction = level > 0 && coarse(halov,0);
	  if (restriction && !is_local(aparent(0)))
	    rcv_pid_append (halo_restriction->snd, aparent(0).pid, point);
	  if (!is_leaf(cell)) {
	    if (!restriction)
	      foreach_neighbor()
		if (is_leaf(cell) && !is_boundary(cell))
		  restriction = true, break;
	    if (restriction) {
	      cell.flags |= halo;
	      halov[] = true;
	    }
	    else {
	      cell.flags &= ~halo;
	      halov[] = false;
	    }
	  }
	}
	else { // non-local cell
	  bool restriction = level > 0 && (aparent(0).flags & halo);
	  if (restriction && is_local(aparent(0)))
	    rcv_pid_append (halo_restriction->rcv, cell.pid, point);
	  if (!is_leaf(cell)) {
	    if (!restriction)
	      // this is necessary so that boundary cells are recognised
	      foreach_neighbor()
		if (allocated(0) && is_leaf(cell) && is_local(cell))
		  restriction = true, break;
	    if (restriction)
	      cell.flags |= halo;
	    else
	      cell.flags &= ~halo;
	  }
	}
	continue; // level == l
      }
      if (is_leaf(cell))
	continue;
    }
    // fixme: optimise, we use only parent values, not neighbors
    mpi_boundary_restriction (mpi_boundary, {halov}, l);
  }
}

trace
void mpi_boundary_refine (scalar * list, int file)
{
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  /* Send refinement cache to each neighboring process. */
  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = quadtree->refined.n;
    MPI_Isend (&quadtree->refined.n, 1, MPI_INT, *dest,
	       REFINE_TAG(), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (quadtree->refined.p, sizeof(Index)/sizeof(int)*len,
		 MPI_INT, *dest, REFINE_TAG(), MPI_COMM_WORLD, &r[nr++]);
  }

  /* Receive refinement cache from each neighboring process. 
   fixme: use non-blocking receives */
  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, REFINE_TAG(),
		    MPI_COMM_WORLD, MPI_STATUS_IGNORE,
		    "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
		      MPI_INT, *source, REFINE_TAG(),
		      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
		      "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
      foreach_cache (refined,)
	if (level <= depth() && allocated(0)) {
	  if (is_leaf(cell)) {
	    bool neighbors = false;
	    foreach_neighbor()
	      if (allocated(0) && (is_active(cell) || is_local(aparent(0))))
		neighbors = true, break;
	    // refine the cell only if it has local neighbors
	    if (neighbors)
	      refine_cell (point, list, 0, &rerefined);
	  }
	}
    }
  }

  /* check that caches were received OK and free ressources */
  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);

  /* update the refinement cache with "re-refined" cells */
  free (quadtree->refined.p);
  quadtree->refined = rerefined;
  
  prof_stop();

  /* if any cell has been re-refined, we repeat the process to take care
     of recursive refinements induced by the 2:1 constraint */
  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list, 0);
}

static void check_depth()
{
#if DEBUG_MPI 
  int max = 0;
  foreach_cell_all()
    if (level > max)
      max = level;
  if (depth() != max) {
    FILE * fp = fopen ("layer", "w");
    fprintf (fp, "depth() = %d, max = %d\n", depth(), max);
    for (int l = 0; l <= depth(); l++) {
      Layer * L = quadtree->L[l];
      fprintf (fp, "Layer level = %d, nc = %d, len = %d\n", l, L->nc, L->len);
      for (int i = 0; i < L->len; i++)
	if (L->m[i]) {
	  fprintf (fp, "  i = %d, refcount = %d\n", i,
		   *((int *)(((char *)L->m[i]) + L->len*sizeof(char *))));
	  for (int j = 0; j < L->len; j++)
	    if (L->m[i][j]) {
	      fprintf (fp, "    j = %d\n", j);
	      fprintf (fp, "point %g %g %d\n",
		       X0 + (i - 1.5)*L0/(1 << l), Y0 + (j - 1.5)*L0/(1 << l),
		       l);
	    }
	}
    }
    fclose (fp);
    fp = fopen ("colls", "w");
    output_cells (fp);
    fclose (fp);
    assert (false);
  }
#endif
}

typedef struct {
  int refined, leaf;
} Remote;

#define REMOTE() ((Remote *)&val(remote,0))

trace
void mpi_boundary_coarsen (int l, int too_fine)
{
  if (npe() == 1)
    return;
  
  check_depth();

  assert (sizeof(Remote) == sizeof(double));
  
  scalar remote[];
  foreach_cell() {
    if (level == l) {
      if (is_local(cell)) {
	REMOTE()->refined = is_refined(cell);
	REMOTE()->leaf = is_leaf(cell);
      }
      else {
	REMOTE()->refined = true;
	REMOTE()->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }
  mpi_boundary_restriction (mpi_boundary, {remote}, l);
  
  foreach_cell() {
    if (level == l) {
      if (!is_local(cell)) {
	if (is_refined(cell) && !REMOTE()->refined)
	  coarsen_cell_recursive (point, NULL);
	else if (is_leaf(cell) && cell.neighbors && REMOTE()->leaf) {
	  int pid = cell.pid;
	  foreach_child()
	    cell.pid = pid;
	}
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }

  check_depth();

  if (l > 0) {
    foreach_cell() {
      if (level == l) {
	remote[] = is_local(cell) ? cell.neighbors : 0;
	continue;
      }
      if (is_leaf(cell))
	continue;
    }
    mpi_boundary_restriction (mpi_boundary, {remote}, l);
    foreach_cell() {
      if (level == l)
	if (!is_local(cell) && is_local(aparent(0)) && remote[]) {
	  aparent(0).flags &= ~too_fine;
	  continue;
	}
      if (is_leaf(cell))
	continue;
    }
  }
}

static void flag_border_cells()
{
  foreach_cell() {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
      foreach_neighbor() {
	if (!is_local(cell) || (level > 0 && !is_local(aparent(0))))
	  flags |= border, break;
	// root cell
	if (is_refined(cell))
	  foreach_child()
	    if (!is_local(cell))
	      flags |= border, break;
	if (flags & border)
	  break;
      }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;
      //      continue; // fixme
    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
	foreach_child()
	  cell.flags &= ~border;
	if (is_border(cell)) {
	  bool remote = false;
	  foreach_neighbor (GHOSTS/2)
	    if (!is_local(cell))
	      remote = true, break;
	  if (remote)
	    foreach_child()
	      cell.flags |= border;
	}
      }
      continue;
    }
  }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}

trace
void mpi_partitioning()
{
  prof_start ("mpi_partitioning");

  long nt = 0;
  foreach()
    nt++;

  /* set the pid of each cell */
  long i = 0;
  quadtree->dirty = true;
  foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
	cell.pid = balanced_pid (i++, nt, npe());
	if (cell.neighbors > 0) {
	  int pid = cell.pid;
	  foreach_child()
	    cell.pid = pid;
	}
	if (!is_local(cell))
	  cell.flags &= ~active;
      }
      else {
	cell.pid = child(0).pid;
	bool inactive = true;
	foreach_child()
	  if (is_active(cell))
	    inactive = false, break;
	if (inactive)
	  cell.flags &= ~active;
      }
    }

  flag_border_cells();
  
  prof_stop();
  
  mpi_boundary_update();
}

/**
# *z_indexing()*: fills *index* with the Z-ordering index.
   
If `leaves` is `true` only leaves are indexed, otherwise all active
cells are indexed. 

On the master process (`pid() == 0`), the function returns the
(global) maximum index (and -1 on all other processes).

On a single processor, we would just need something
like (for leaves)

~~~literatec
double i = 0;
foreach()
  index[] = i++;
~~~

In parallel, this is a bit more difficult. */

trace
double z_indexing (scalar index, bool leaves)
{
  /**
  ## Size of subtrees

  We first compute the size of each subtree. We use *index* to store
  this. */
  
  scalar size[];

  /**
  The size of leaf "subtrees" is one. */

  foreach()
    size[] = 1;
  
  /**
  We do a (parallel) restriction to compute the size of non-leaf
  subtrees. */

  // fixme: mpi_boundary_restriction should be enough
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
  fixme: doc */
  
  foreach_level(0)
    index[] = 0;
  for (int l = 0; l < depth(); l++) {
    // fixme: mpi_boundary_restriction should be enough
    boundary_iterate (restriction, {index}, l);
    foreach_cell() {
      if (level == l) {
	if (is_leaf(cell)) {
	  if (is_local(cell) && cell.neighbors) {
	    int i = index[];
	    foreach_child()
	      index[] = i;
	  }
	}
	else { // not leaf
	  bool loc = is_local(cell);
	  if (!loc)
	    foreach_child()
	      if (is_local(cell))
		loc = true, break;
	  if (loc) {
	    int i = index[] + !leaves;
	    foreach_child() {
	      index[] = i;
	      i += size[]; 
	    }
	  }
	}
	continue;
      }
      if (is_leaf(cell))
	continue;
    }
  }
  // fixme: mpi_boundary_restriction should be enough
  boundary_iterate (restriction, {index}, depth());

  return maxi;
}
