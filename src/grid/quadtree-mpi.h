// fixme: this does not work with dynamic quadtrees

#define DEBUG 0

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
      foreach_cache_level(m->snd_halo[i][l], l,)
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
      foreach_cache_level(m->rcv_halo[i][l], l,) {
#if DEBUG
	fprintf (stderr, "%g %g rcv\n", x, y);
#endif
	for (scalar s in list)
	  s[] = *b++;
      }
    }

  /* check that ghost values where received OK and free send buffers */
  for (int i = 0; i < m->nsnd; i++)
    if (m->snd_halo[i][l].n > 0) {
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
  add_boundary (mpi_boundary);
}

#define DEBUG 1

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
  scalar pid[];
  int i = 0;
  foreach_cell_post (is_active (cell) || cell.neighbors > 0)
    if (is_active (cell)) {
      if (is_leaf (cell))
	pid[] = i++/nf;
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
  int adjacency_rcv[size];
  for (int i = 0; i < size; i++)
    adjacency_rcv[i] = 0;

  /* we build arrays of ghost cell indices for restriction */
  ((Quadtree *)grid)->dirty = true;
  foreach_cell() {
    if (!is_active(cell)) {
      if (!is_leaf(cell)) {
	// decrement_neighbors (&point);
	for (int k = -GHOSTS/2; k <= GHOSTS/2; k++)
	  for (int l = -GHOSTS/2; l <= GHOSTS/2; l++)
	    if (!is_active(neighbor(k,l))) {
	      neighbor(k,l).neighbors--;
	      if (neighbor(k,l).neighbors == 0)
		free_children((Quadtree *)grid,k,l);
	    }
      }
      int nactive = 0;
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  if (is_active(neighbor(o,p)))
	    nactive++;
      if (nactive > 0) {
	int i;
	for (i = 0; i < m->nrcv; i++)
	  if (pid[] == m->rcv[i])
	    break;
	if (i == m->nrcv) {
	  m->rcv = realloc (m->rcv, ++m->nrcv*sizeof (int));
	  m->rcv[m->nrcv-1] = pid[];
	  m->rcv_halo = realloc (m->rcv_halo, m->nrcv*sizeof (CacheLevel *));
	  m->rcv_halo[m->nrcv-1] = malloc ((depth() + 1)*sizeof (CacheLevel));
	  for (int j = 0; j <= depth(); j++)
	    cache_level_init (&m->rcv_halo[m->nrcv-1][j]);
	}
	cache_level_append (&m->rcv_halo[i][level], point);
	adjacency_rcv[(int)pid[]]++;
      }
    }
    if (is_leaf(cell))
      continue;
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
  void output_cells (FILE * fp);

  char name[80];
  sprintf (name, "halo-%d", pid());
  FILE * fp = fopen (name, "w");

  // local halo
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l)
      foreach_child()
        fprintf (fp, "%g %g %g %d\n", x, y, pid[], level);
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
      fprintf (fp, "%g %g %g %d %d\n", x, y, pid[], level, cell.neighbors);
  fclose (fp);

  sprintf (name, "ghost-%d", pid());
  fp = fopen (name, "w");
  for (int i = 0; i < m->nsnd; i++)
    for (int j = 0; j <= depth(); j++)
      foreach_cache_level (m->snd_halo[i][j], j,)
	fprintf (fp, "%g %g %d %d\n", x, y, m->snd[i], cell.neighbors);
  fclose (fp);
#endif
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

