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
#if 0
  if (is_refined(cell))
    level++;
#endif
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

static Rcv * rcv_pid_pointer (RcvPid * p, Array * a, int pid)
{
  assert (pid >= 0 && pid < npe());
  int * apid = (int *) a->p;
  for (int i = 0; i < a->len/sizeof(int); i++)
    if (apid[i] == pid)
      return NULL;
  array_append (a, &pid, sizeof(int));
  
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

static void rcv_pid_append (RcvPid * p, Array * a, int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, a, pid);
  if (rcv)
    rcv_append (point, rcv);
}

static Rcv * rcv_pid_pointer1 (RcvPid * p, int pid)
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

static const unsigned short used = 1 << user;

static void rcv_pid_append1 (RcvPid * p, int pid, Point point)
{
  Rcv * rcv = rcv_pid_pointer1 (p, pid);
  if (rcv) {
    cell.flags |= used;
    rcv_append (point, rcv);
  }
}

#if 0
static void rcv_pid_append_children (RcvPid * p, Array * a, int pid,
				     Point point)
{
  Rcv * rcv = rcv_pid_pointer (p, a, pid);
  if (rcv)
    foreach_child()
      rcv_append (point, rcv);
}
#endif

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

static size_t apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
			bool children)
{
  double * b = rcv->buf;
  if (!children)
    foreach_cache_level(rcv->halo[l], l,) {
      for (scalar s in list)
	s[] = *b++;
      for (vector v in listf) {
	foreach_dimension() {
	  if (allocated(-1))
	    v.x[] = *b;
	  b++;
	  if (allocated(1))
	    v.x[1] = *b;
	  b++;
	}
      }
    }
  else
    foreach_cache_level(rcv->halo[l], l,)
      if (is_refined(cell))
	foreach_child()
	  for (scalar s in list)
	    s[] = *b++;
  size_t size = b - (double *) rcv->buf;
  free (rcv->buf);
  rcv->buf = NULL;
  return size;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l,
			     bool children)
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
      long n = 0;
      if (!children)
	n = rcv->halo[l].n;
      else
	foreach_cache_level(rcv->halo[l], l,)
	  if (is_refined(cell))
	    n += (1 << dimension);
      if (n > 0) {
	rcv->buf = malloc (sizeof (double)*n*len);
#if 0
	fprintf (stderr, "%s receiving %d doubles from %d level %d\n",
		 m->name, n*len, rcv->pid, l);
	fflush (stderr);
#endif
#if 1 /* initiate non-blocking receive */
	MPI_Irecv (rcv->buf, n*len, MPI_DOUBLE, rcv->pid,
		   BOUNDARY_TAG(l), MPI_COMM_WORLD, &r[nr]);
	rrcv[nr++] = rcv;
#else /* blocking receive (useful for debugging) */
	MPI_Recv (rcv->buf, n*len, MPI_DOUBLE, rcv->pid,
		  BOUNDARY_TAG(l), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	apply_bc (rcv, list, listf, l, children);
#endif
      }
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
      size_t size = apply_bc (rcv, list, listf, l, children);
      int rlen;
      MPI_Get_count (&s, MPI_DOUBLE, &rlen);
      assert (rlen == size);
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

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l,
			  bool children)
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
      double * b;
      long n = 0;
      if (!children) {
	n = rcv->halo[l].n;
	rcv->buf = malloc (sizeof (double)*n*len);
	b = rcv->buf;
	foreach_cache_level(rcv->halo[l], l,) {
	  for (scalar s in list)
	    *b++ = s[];
	  for (vector v in listf)
	    foreach_dimension() {
	      *b++ = v.x[];
	      *b++ = allocated(1) ? v.x[1] : undefined;
	    }
	}
      }
      else { // send children of refined cells
	foreach_cache_level(rcv->halo[l], l,)
	  if (is_refined(cell))
	    n += (1 << dimension);
	if (n > 0) {
	  rcv->buf = malloc (sizeof(double)*n*len);
	  b = rcv->buf;
	  foreach_cache_level(rcv->halo[l], l,)
	    if (is_refined(cell))
	      foreach_child()
		for (scalar s in list)
		  *b++ = s[];
	}
      }
#if 0
      fprintf (stderr, "%s sending %d doubles to %d level %d\n",
	       m->name, rcv->halo[l].n*len, rcv->pid, l);
      fflush (stderr);
#endif
      if (n > 0)
	MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
		   MPI_DOUBLE, rcv->pid, BOUNDARY_TAG(l), MPI_COMM_WORLD,
		   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l, bool children)
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
  rcv_pid_send (m->snd, listr, listf, l, children);
  rcv_pid_receive (m->rcv, listr, listf, l, children);
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
  free (m);
}

trace
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l, false);
  rcv_pid_sync (&m->restriction_root, list, l, false);
}

trace
void mpi_boundary_restriction_children (const Boundary * b, scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l, true);
  rcv_pid_sync (&m->restriction_root, list, l, true);
}

trace
static void mpi_boundary_halo_restriction (const Boundary * b,
					   scalar * list, int l)
{
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->halo_restriction, list, l, false);
}

trace
static void mpi_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, int l, int depth)
{
  MpiBoundary * m = (MpiBoundary *) b;
  if (true /*l == depth*/) {
    rcv_pid_sync (&m->restriction, list, l, false);
    rcv_pid_sync (&m->restriction_root, list, l, false);
  }
  else
    rcv_pid_sync (&m->prolongation, list, l, false);
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
    sprintf (name, "mpi-rborder-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
    sprintf (name, "coarsened-%d", pid()); remove (name);
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
	       prefix, x, y, z, level, cell.neighbors,
	       (cell.flags & ignored) ? npe() : cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-rborder", prefix);
  foreach_cell()
    if (is_rborder(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
	       prefix, x, y, z, level, cell.neighbors,
	       (cell.flags & ignored) ? npe() : cell.pid);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
  foreach_cell() {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
	       prefix, x, y, z, level, cell.neighbors,
	       (cell.flags & ignored) ? npe() : cell.pid, cell.flags & leaf);
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

  fp = fopen_prefix (fp1, "coarsened", prefix);
  for (int l = 0; l < depth(); l++)
    foreach_cache_level (quadtree->coarsened, l,)
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
  int i = p.i - point.i, j = p.j - point.j;
  if (j == 0) {
    for (int j = -1; j <= 1; j += 2)
      if (is_refined(neighbor(0,j)))
	return true;
  }
  else if (i == 0) {
    for (int i = -1; i <= 1; i += 2)
      if (is_refined(neighbor(i,0)))
	return true;
  }
  else if (is_refined(neighbor(i,0)) || is_refined(neighbor(0,j)))
    return true;
  return false;
}
		
static int locals_pids (Point point, Array * pids)
{
  if (is_leaf(cell)) { // prolongation
    if (is_local(cell)) {
      Point p = point;
      foreach_neighbor(1) {
	if (cell.pid >= 0 && !is_local(cell) &&
	    (is_refined(cell) || is_local_prolongation (point, p)))
	  append_pid (pids, point);
	if (is_refined(cell))
	  foreach_child()
	    if (cell.pid >= 0 && !is_local(cell))
	      append_pid (pids, point);
      }
    }
  }
  else
    foreach_neighbor(1) {
      if (cell.pid >= 0 && !is_local(cell))
	append_pid (pids, point);
      if (is_refined(cell))
	foreach_child()
	  if (cell.pid >= 0 && !is_local(cell))
	    append_pid (pids, point);
    }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{
  foreach_child()
    if (cell.pid >= 0 && !is_local(cell))
      append_pid (pids, point);
  return pids->len/sizeof(int);
}

trace
void mpi_boundary_update_new()
{
  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * restriction = &m->restriction;
  SndRcv * restriction_root = &m->restriction_root;

#if 0
  char name[80];
  sprintf (name, "update-%d", pid());
  FILE * fp = fopen (name, "w");
#endif
  
  snd_rcv_free (restriction);
  snd_rcv_free (restriction_root);

  foreach_cell() {
    if (cell.neighbors) {
      // sending
      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
	foreach_child()
	  if (is_local(cell))
	    for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
	      rcv_pid_append1 (restriction->snd, *p, point);
	free (pids.p);
      }
      // receiving
      bool locals = false;
      if (is_leaf(cell)) { // prolongation
	if (cell.pid >= 0 && !is_local(cell)) {
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
	  if (cell.pid >= 0 && !is_local(cell))
            rcv_pid_append1 (restriction->rcv, cell.pid, point);

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
		  rcv_pid_append1 (restriction_root->snd, *p, point);
	    free (pids.p);
	  }
	}
	// receiving
	else if (cell.pid >= 0) {
	  bool root = false;
	  foreach_child()
	    if (is_local(cell))
	      root = true, break;
	  if (root) {
	    int pid = cell.pid;
	    foreach_neighbor()
	      if (cell.pid >= 0 && !is_local(cell))
		rcv_pid_append1 (restriction_root->rcv, pid, point);
	  }
	}
      }
    }
#if 0
    // we mark all non-local cells as "ignored"
    if (cell.pid >= 0 && !is_local(cell))
      cell.flags |= ignored;
    else // local cells are not ignored
      cell.flags &= ~ignored;
#endif
  }

  // all received non-local cells are marked as "not ignored"
  
  
  /* we mark/remove unused cells
     we do a breadth-first traversal from fine to coarse, so that
     coarsening of unused cells can proceed fully. See also
     quadtree-common.h:coarsen_function() */
  for (int l = depth(); l >= 0; l--)
    foreach_cell()
      if (level == l) {
	if (cell.flags & used)
	  cell.flags &= ~used;
	else if (!is_local(cell) && !is_prolongation(cell)) {
	  cell.pid = npe();
	  if (is_leaf(cell) && cell.neighbors)
	    foreach_child()
	      cell.pid = npe();
	}
	if (is_refined(cell)) {
	  bool coarsen = true;
	  foreach_child()
	    if (!is_leaf(cell) || cell.pid != npe())
	      coarsen = false, break;
	  if (coarsen)
	    coarsen_cell (point, NULL, NULL);
	}
	continue; // level == l
      }
  
#if 0
  fclose (fp);
#endif
  
}
  
trace
void mpi_boundary_update()
{
  if (npe() == 1)
    return;

  prof_start ("mpi_boundary_update");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;  

#if 0
  SndRcv * restriction = &m->restriction;
  SndRcv * prolongation = &m->prolongation;

  snd_rcv_free (restriction);
  snd_rcv_free (prolongation);
  
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
	  Array apro = {NULL, 0, 0}, ares = {NULL, 0, 0};
	  Point p = point;
	  foreach_neighbor()
	    if (cell.pid >= 0) {
	      if (!is_local(cell)) {
		rcv_pid_append (res, &ares, cell.pid, p);
		if (is_leaf(cell) || is_prolongation(cell))
		  rcv_pid_append (pro, &apro, cell.pid, p);
	      }
	      if (level > 0 && !is_local(aparent(0)))
		rcv_pid_append (res, &ares, aparent(0).pid, p);
	      // root cell: fixme: optimise?
	      if (is_refined(cell))
		foreach_child()
		  if (!is_local(cell))
		    rcv_pid_append (res, &ares, cell.pid, p);
	    }
	  if (is_leaf(cell)) {
	    foreach_neighbor (1)
	      if (cell.pid >= 0 && !is_local(cell))
		rcv_pid_append (pro, &apro, cell.pid, p);
	    if (cell.neighbors) {
	      // prolongation
	      Array cpro = {NULL, 0, 0}, cres = {NULL, 0, 0};
	      foreach_neighbor (1)
		if (cell.pid >= 0 && cell.neighbors)
		  foreach_child()
		    if (!is_local(cell)) {
		      rcv_pid_append (res, &ares, cell.pid, p);
		      rcv_pid_append_children (res, &cres, cell.pid, p);
		      if (!is_refined(cell)) {
			rcv_pid_append (pro, &apro, cell.pid, p);
			rcv_pid_append_children (pro, &cpro, cell.pid, p);
		      }
		    }
	      free (cpro.p); free (cres.p);
	    }
	  }
	  cell.flags &= ~ignored;
	  free (apro.p); free (ares.p);
	}
	else {
	  // =========================================
	  // non-local cell: do we need to receive it?
	  RcvPid * pro = prolongation->rcv;
	  RcvPid * res = restriction->rcv;
	  bool used = false, coarsen = true;
	  Array apro = {NULL, 0, 0}, ares = {NULL, 0, 0};
	  int pid = cell.pid;
	  Point p = point;
	  foreach_neighbor()
	    if (allocated(0) && !is_boundary(cell)) {
	      if (is_local(cell)) {
		if (!used)
		  used = true, rcv_pid_append (res, &ares, pid, p);
		if (is_leaf(cell) || !is_active(cell))
		  rcv_pid_append (pro, &apro, pid, p),
		    break;
	      }
	      else if (!used) {
		// root cell: fixme: optimise?
		if (is_refined(cell))
		  foreach_child()
		    if (is_local(cell))
		      used = true, rcv_pid_append (res, &ares, pid, p), break;
		/* see the corresponding condition in
		   mpi_boundary_refine(), and src/test/mpi-refine1.c
		   on why `level > 0 && is_local(aparent(0)))` is necessary. */
		if (level > 0 && is_local(aparent(0)))
		  used = true, rcv_pid_append (res, &ares, pid, p);
	      }
	    }
	  if (is_leaf(cell)) {
	    foreach_neighbor (1)
	      if (allocated(0) && is_local(cell))
		rcv_pid_append (pro, &apro, pid, p),
		  break;
	    if (cell.neighbors) {
	      // prolongation
	      Array cpro = {NULL, 0, 0}, cres = {NULL, 0, 0};
	      bool used = false;
	      foreach_neighbor (1)
		if (allocated(0) && is_active(cell) && cell.neighbors)
		  foreach_child()
		    if (is_local(cell)) {
		      used = true;
		      rcv_pid_append (res, &ares, pid, p);
		      rcv_pid_append_children (res, &cres, pid, p);
		      if (!is_refined(cell)) {
			rcv_pid_append (pro, &apro, pid, p);
			rcv_pid_append_children (pro, &cpro, pid, p);
			break;
		      }
		    }
	      free (cpro.p); free (cres.p);
	      if (used)
		foreach_child()
		  cell.flags |= rborder;
	      else
		foreach_child()
		  cell.flags &= ~rborder;
	    }
	  }
	  else { // not leaf
	    if (used)
	      foreach_child()
		cell.flags &= ~ignored;
	    else {
	      if (coarsen)
		foreach_child()
		  if (!(cell.flags & ignored))
		    coarsen = false, break;
	      // non-local coarse cell with no local/active neighbors:
	      // destroy its children
	      if (coarsen)
		coarsen_cell (point, NULL, NULL);
	    }
	  } // not leaf
	  if (used) {
	    cell.flags &= ~ignored;
	    if (is_leaf(cell) && cell.neighbors)
	      foreach_child()
		cell.flags &= ~ignored;
	    cell.flags |= rborder;
	  }
	  else {
	    cell.flags |= ignored;
	    if (is_leaf(cell) && cell.neighbors)
	      foreach_child()
		cell.flags |= ignored;
	    cell.flags &= ~rborder;
	  }
	  free (apro.p); free (ares.p);
	}
	continue; // level == l
      }
    
      if (is_leaf(cell))
	continue;
    }

#endif  
  
  mpi_boundary_update_new();
  
  prof_stop();  

#if DEBUG_MPI
  debug_mpi (NULL);
#endif

#if 1
  // halo restriction
  SndRcv * halo_restriction = &m->halo_restriction;
  snd_rcv_free (halo_restriction);
  scalar halov[];
  quadtree->dirty = true;
  for (int l = 0; l <= depth(); l++) {
    Array apro = {NULL, 0, 0};
    foreach_cell() {
      if (level == l) {
	if (is_local(cell)) {
	  bool restriction = level > 0 && coarse(halov,0);
	  if (restriction && !is_local(aparent(0))) {
	    RcvPid * halo_res = halo_restriction->snd;
	    apro.len = 0;
	    rcv_pid_append (halo_res, &apro, aparent(0).pid, point);
	  }
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
	  if (restriction && is_local(aparent(0))) {
	    RcvPid * halo_res = halo_restriction->rcv;
	    apro.len = 0;
	    rcv_pid_append (halo_res, &apro, cell.pid, point);
	  }
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
    free (apro.p);
    // fixme: optimise, we use only parent values, not neighbors
    mpi_boundary_restriction (mpi_boundary, {halov}, l);
  }
#endif  
}

trace
void mpi_boundary_refine (scalar * list)
{
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  char name[80];
  sprintf (name, "refine-%d", pid());
  FILE * fp = fopen (name, "w");
  
  /* Send refinement cache to each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  MPI_Request r[2*snd->npid];
  int nr = 0;
  for (int i = 0; i < snd->npid; i++) {
    int pid = snd->rcv[i].pid, len = quadtree->refined.n;
    MPI_Isend (&quadtree->refined.n, 1, MPI_INT, pid,
	       REFINE_TAG(), MPI_COMM_WORLD, &r[nr++]);
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
	      if (allocated(0) && (is_active(cell) || is_local(aparent(0))))
		neighbors = true, break;
	    // refine the cell only if it has local neighbors
	    fprintf (fp, "%g %g\n", x, y);
	    if (neighbors) {
	      refine_cell (point, list, 0, &rerefined, NULL);
	    }
	  }
	}
    }
  }

  fclose (fp);
  
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
    mpi_boundary_refine (list);
}

static void check_depth()
{
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
}

trace
void mpi_boundary_coarsen (int l, int too_fine)
{
  //  prof_start ("mpi_boundary_coarsen");

  char name[80];
  sprintf (name, "coarsen-%d-%d", pid(), l);
  FILE * fp = fopen (name, "w");

#if 1
  check_depth();
  
  scalar refined[];
  quadtree_trash ({refined}); // fixme
  foreach_cell() {
    if (level == l) {
      if (is_local(cell))
	refined[] = is_refined(cell);
      continue;
    }
    if (is_leaf(cell))
      continue;
  }
  mpi_boundary_restriction (mpi_boundary, {refined}, l);

  foreach_cell() {
    if (level == l) {
      if (!is_local(cell) && is_refined(cell) &&
	  !isnan(refined[]) && !refined[]) {
	fprintf (fp, "%g %g %d refined\n", x, y, cell.pid);
	//	coarsen_cell_recursive (point, NULL, NULL);
	assert (coarsen_cell (point, NULL, NULL));
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }

  check_depth();

  if (l > 0) {
    // fixme: optimize cell.neighbors
    foreach_cell() {
      if (level == l) {
	if (is_local(cell))
	  refined[] = cell.neighbors;
	continue;
      }
      if (is_leaf(cell))
	continue;
    }
    mpi_boundary_restriction (mpi_boundary, {refined}, l);
    foreach_cell() {
      if (level == l) {
	if (!is_local(cell) && !isnan(refined[]) && refined[])
	  aparent(0).flags &= ~too_fine;
	continue;
      }
      if (is_leaf(cell))
	continue;
    }
  }
#else
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  /* Send coarsening cache to each neighboring process. */
  RcvPid * snd = mpi->restriction.snd;
  MPI_Request r[2*snd->npid];
  int nr = 0;
  for (int i = 0; i < snd->npid; i++)
    if (snd->rcv[i].maxdepth > l) {
      int pid = snd->rcv[i].pid, len = quadtree->coarsened.n;
      MPI_Isend (&quadtree->coarsened.n, 1, MPI_INT, pid, COARSEN_TAG(l),
		 MPI_COMM_WORLD, &r[nr++]);
      fprintf (fp, "%d sending %d to %d\n", pid(), len, pid);
      if (len > 0) {
	MPI_Isend (quadtree->coarsened.p, sizeof(IndexLevel)/sizeof(int)*len,
		   MPI_INT, pid, COARSEN_TAG(l), MPI_COMM_WORLD, &r[nr++]);
      }
    }
  
  /* Receive coarsening cache from each (fine enough) neighboring process. 
   fixme: use non-blocking receives */
  RcvPid * rcv = mpi->restriction.rcv;
  for (int i = 0; i < rcv->npid; i++)
    if (rcv->rcv[i].maxdepth > l) {
      int pid = rcv->rcv[i].pid, len;
      MPI_Recv (&len, 1, MPI_INT, pid, COARSEN_TAG(l),
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      fprintf (fp, "%d receiving %d from %d\n", pid(), len, pid);
      if (len > 0) {
	IndexLevel p[len];
	MPI_Recv (p, sizeof(IndexLevel)/sizeof(int)*len,
		  MPI_INT, pid, COARSEN_TAG(l),
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	CacheLevel coarsened = {p, len, len};
	foreach_cache_level (coarsened, l,)
	  if (allocated(0) && is_refined(cell)) {
	    fprintf (fp, "%g %g %d coarsened\n", x, y, pid);
	    coarsen_cell_recursive (point, NULL, NULL);
	  }
      }
    }

  /* check that caches were received OK and free resources */
  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);
#endif
  
  fclose (fp);
  
  //  prof_stop();
}

static void flag_border_cells()
{
  foreach_cell() {
    if (is_active(cell)) {
      short flags = cell.flags & ~(border|rborder);
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
      cell.flags &= ~(border|rborder);
      //      continue; // fixme
    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
	foreach_child()
	  cell.flags &= ~(border|rborder);
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

static int balanced_pid (long index, long nt)
{  
  long ne = max(1, nt/npe()), nr = nt % npe();
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(npe() - 1, pid);
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
	cell.pid = balanced_pid (i++, nt);
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
  boundary_iterate (restriction, {index}, depth());

  return maxi;
}
