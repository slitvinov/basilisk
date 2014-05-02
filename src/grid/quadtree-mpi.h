void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit ((void (*)(void)) MPI_Finalize);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      stdout = freopen (name, "w", stdout);
      sprintf (name, "log-%d", mpi_rank);
      stderr = freopen (name, "w", stderr);
    }
  }
}

typedef struct {
  Boundary parent;
  
  int nrcv;   // the number of PEs to receive ghost values from
  int * rcv;  // the rank of each PE
  // ghost cell indices for each PE and each level
  CacheLevel ** rcv_halo;

  int nsnd;   // the number of PEs to send ghost values to
  int * snd;  // the rank of each PE
  // ghost cell indices for each PE and each level
  CacheLevel ** snd_halo;
} MpiBoundary;

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  for (int i = 0; i < m->nrcv; i++)
    free_cache (m->rcv_halo[i]);
  free (m->rcv_halo);
  free (m->rcv);
  for (int i = 0; i < m->nsnd; i++)
    free_cache (m->snd_halo[i]);
  free (m->snd_halo);
  free (m->snd);
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
  for (int i = 0; i < m->nsnd; i++) 
    if (m->snd_halo[i][l].n > 0) {
      buf[i] = malloc (sizeof (double)*m->snd_halo[i][l].n*len);
      double * b = buf[i];
      foreach_cache_level(m->snd_halo[i][l], l)
	for (scalar s in list)
	  *b++ = s[];
      MPI_Isend (buf[i], m->snd_halo[i][l].n*len, MPI_DOUBLE, m->snd[i], l, 
		 MPI_COMM_WORLD, &r[i]);
    }

  /* receive ghost values */
  for (int i = 0; i < m->nrcv; i++) 
    if (m->rcv_halo[i][l].n > 0) {
      double buf[m->rcv_halo[i][l].n*len], * b = buf;
      MPI_Status s;
      MPI_Recv (buf, m->rcv_halo[i][l].n*len, MPI_DOUBLE, m->rcv[i], l, 
		MPI_COMM_WORLD, &s);
      foreach_cache_level(m->rcv_halo[i][l], l)
	for (scalar s in list)
	  s[] = *b++;
    }

  /* check that ghost values where received OK and free send buffers */
  for (int i = 0; i < m->nsnd; i++)
    if (m->snd_halo[i][l].n > 0) {
      MPI_Status s;
      MPI_Wait (&r[i], &s);
      free (buf[i]);
    }
}

static void mpi_boundary_halo_prolongation (const Boundary * b,
					    scalar * list, int l)
{}

static Boundary * mpi_boundary = NULL;

void mpi_boundary_new()
{
  mpi_boundary = calloc (1, sizeof (MpiBoundary));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->halo_restriction = mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation = mpi_boundary_halo_prolongation;
  add_boundary (mpi_boundary);
}

#define DEBUG 0

void mpi_partitioning()
{
  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  int nf = 0;
  foreach()
    nf++;
  
  int size;
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  nf = nf/size + 1;

  // set the pid of each cell
  int i = 0;
  foreach_cell_post (is_active (cell) || cell.neighbors > 0) {
    if (is_active (cell)) {
      if (is_leaf (cell))
	cell.pid = i++/nf;
      else {
	cell.pid = -1;
	for (int i = 0; i <= 1; i++)
	  for (int j = 0; j <= 1; j++)
	    if (cell.pid == -1)
	      cell.pid = child(i,j).pid;
	    else if (child(i,j).pid == pid())
	      cell.pid = pid();
      }
    }
    cell.neighbors = 0;
  }

  foreach_halo_coarse_to_fine(depth())
    cell.pid = aparent(0,0).pid;

  // set the number of active neighboring cells belonging to the
  // current process
  foreach_cell() {
    if (cell.pid != pid() || !is_active (cell))
      continue;
    else {
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  neighbor(o,p).neighbors++;
    }
  }

  // Remove cells which do not have any neighbors
  foreach_cell() {
    if (cell.neighbors == 0) {
      free (allocated(0,0));
      allocated(0,0) = NULL;
      point.back->dirty = true;
    }
    else
      cell.neighbors = 0;
  }

  // set the number of leaf neighboring cells
  foreach_cell()
    if (is_leaf (cell)) {
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  if (allocated(o,p))
	    neighbor(o,p).neighbors++;
      continue;
    }

#if DEBUG
  char name[80];
  sprintf (name, "log-%d", pid());
  FILE * logfile = fopen (name, "w");

  sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  // local halo
  foreach_halo_coarse_to_fine(depth())
    fprintf (logfile, "%g %g %d %d\n", x, y, cell.pid, level);

  // local restriction
  sprintf (name, "restrict-%d", pid());
  fp = fopen (name, "w");
  foreach_halo_fine_to_coarse()
    fprintf (fp, "%g %g %d %d\n", x, y, cell.pid, level);
  fclose (fp);

  // MPI communication for leaf-level restriction
  sprintf (name, "mpi-%d", pid());
  fp = fopen (name, "w");  
  foreach_cell() {
    if (!is_active (cell))
      continue;
    else if (cell.pid != pid() && cell.neighbors > 0)
      fprintf (fp, "%g %g %d %d\n", x, y, cell.pid, level);
  }
  fclose (fp);
#endif

  /* this is the adjacency vector i.e. adjacency_rcv[x] is different from
     zero if we need to receive data from proc x */
  int adjacency_rcv[size];
  for (int i = 0; i < size; i++)
    adjacency_rcv[i] = 0;

  /* we build arrays of ghost cell indices for leaf-level restriction */
  foreach_cell() {
    if (!is_active (cell))
      continue;
    else if (cell.pid != pid() && cell.neighbors > 0) {
      int i;
      for (i = 0; i < m->nrcv; i++)
	if (cell.pid == m->rcv[i])
	  break;
      if (i == m->nrcv) {
	m->rcv = realloc (m->rcv, ++m->nrcv*sizeof (int));
	m->rcv[m->nrcv-1] = cell.pid;
	m->rcv_halo = realloc (m->rcv_halo, m->nrcv*sizeof (CacheLevel *));
	m->rcv_halo[m->nrcv-1] = malloc ((depth() + 1)*sizeof (CacheLevel));
	for (int j = 0; j <= depth(); j++)
	  cache_level_init (&m->rcv_halo[m->nrcv-1][j]);
      }
      cache_level_append (&m->rcv_halo[i][level], point);
      adjacency_rcv[cell.pid]++;
    }
  }

  /* we do a global transpose of the adjacency vector
     i.e. adjacency_snd[x] is different from zero if we need to send
     data to proc x */
  int adjacency_snd[size];
  MPI_Alltoall (adjacency_rcv, 1, MPI_INT, adjacency_snd, 1, MPI_INT,
		MPI_COMM_WORLD);

  /* we send the ghost cell indices to their respective PEs */
  for (int i = 0; i < m->nrcv; i++)
    for (int j = 0; j <= depth(); j++) {
      CacheLevel halo = m->rcv_halo[i][j];
      MPI_Request r;
      MPI_Isend (&halo.n, 1, MPI_INT, m->rcv[i], j, MPI_COMM_WORLD, &r);
      if (halo.n > 0)
	MPI_Isend (halo.p, 2*halo.n, MPI_INT, m->rcv[i], j, MPI_COMM_WORLD, &r);
    }

  /* we receive the ghost cell indices from their respective PEs */
  for (int i = 0; i < size; i++)
    if (adjacency_snd[i]) {
      m->snd = realloc (m->snd, ++m->nsnd*sizeof (int));
      m->snd[m->nsnd-1] = i;
      m->snd_halo = realloc (m->snd_halo, m->nsnd*sizeof (CacheLevel *));
      m->snd_halo[m->nsnd-1] = malloc ((depth() + 1)*sizeof (CacheLevel));
      for (int j = 0; j <= depth(); j++) {
	CacheLevel * halo = &m->snd_halo[m->nsnd-1][j];
	MPI_Status s;
	MPI_Recv (&halo->n, 1, MPI_INT, i, j, MPI_COMM_WORLD, &s);
	if (halo->n > 0) {
	  halo->p = malloc (halo->n*sizeof(IndexLevel));
	  MPI_Recv (halo->p, 2*halo->n, MPI_INT, i, j, MPI_COMM_WORLD, &s);
	}
      }
    }

#if DEBUG
  sprintf (name, "ghost-%d", pid());
  fp = fopen (name, "w");
  for (int i = 0; i < m->nsnd; i++)
    for (int j = 0; j <= depth(); j++)
      foreach_cache_level (m->snd_halo[i][j], j)
	fprintf (fp, "%g %g %d %d\n", x, y, m->snd[i], cell.neighbors);
  fclose (fp);
#endif
}
