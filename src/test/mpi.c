@def foreach_cache(_cache,_l) {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  OMP_PARALLEL()
  Quadtree point = *((Quadtree *)grid); point.back = grid;
  point.level = _l;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
    point.j = _cache.p[_k].j;
    POINT_VARIABLES;
@
@define end_foreach_cache() } OMP_END_PARALLEL() }

int * rcv = NULL, nrcv = 0;
CacheLevel ** rcv_halo = NULL;

int * snd = NULL, nsnd = 0;
CacheLevel ** snd_halo = NULL;  

void boundary_mpi (scalar a, int l)
{
  /* send ghost values for a */
  for (int i = 0; i < nsnd; i++) 
    if (snd_halo[i][l].n > 0) {
      double buf[snd_halo[i][l].n];
      int k = 0;
      foreach_cache(snd_halo[i][l], l)
	buf[k++] = a[];
      MPI_Request r;
      MPI_Isend (buf, snd_halo[i][l].n, MPI_DOUBLE, snd[i], l, 
		 MPI_COMM_WORLD, &r);
    }

  /* receive ghost values for a */
  for (int i = 0; i < nrcv; i++) 
    if (rcv_halo[i][l].n > 0) {
      double buf[rcv_halo[i][l].n];
      MPI_Status s;
      MPI_Recv (buf, rcv_halo[i][l].n, MPI_DOUBLE, rcv[i], l, 
		MPI_COMM_WORLD, &s);
      int k = 0;
      foreach_cache(rcv_halo[i][l], l)
	a[] = buf[k++];
    }
}

void halo_restriction_mpi (int l, scalar s)
{
  char name[80];
  sprintf (name, "res-%d", pid());
  FILE * fp = fopen (name, "w");
  if (l < 0) l = depth();
  boundary_mpi (s, l);
  for (l--; l >= 0; l--) {
    CacheLevel halo = ((Quadtree *)grid)->halo[l];
    foreach_cache (halo, l)
      if (is_active (cell)) {
	fprintf (fp, "%g %g %d\n", x, y, level);
	s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
      }
    boundary_mpi (s, l);
  }
  fclose (fp);
}

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  x -= 0.1; y -= 0.1;
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

int main(int argc, char * argv[])
{
  MPI_Init (&argc, &argv);
  MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

  int depth = argc > 1 ? atoi(argv[1]) : 4;
  X0 = Y0 = -0.5;
  init_grid(1);
  while (refine_function (refine_circle, &depth, NULL));

  int nf = 0;
  foreach()
    nf++;
  
  int size;
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  nf = nf/size + 1;

  int i = 0, rank = pid();
  foreach_cell_post (is_active (cell)) {

    /* We set the pid of each cell. */

    if (is_leaf (cell))
      cell.pid = i++/nf;
    else {
      cell.pid = -1;
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
	  if (cell.pid == -1)
	    cell.pid = child(i,j).pid;
	  else if (child(i,j).pid == rank)
	    cell.pid = rank;
    }
  }
    
  foreach_halo()
    cell.pid = aparent(0,0).pid;

  foreach_cell_post (is_active (cell)) {
    /* each cell which does not belong to the PE is set to inactive
       and its neighbors are informed of its change of state. */

    if (cell.pid != rank) {
#if 1
      /* update neighborhood */
      if (is_leaf (cell))
	for (int o = -GHOSTS; o <= GHOSTS; o++)
	  for (int p = -GHOSTS; p <= GHOSTS; p++)
	    neighbor(o,p).neighbors--;
#endif
      cell.flags &= ~(active|leaf);
      point.back->dirty = true;
    }
  }

  /* this is the adjacency vector i.e. adjacency_rcv[x] is different from
     zero if we need to receive data from proc x */
  int adjacency_rcv[size];
  for (int i = 0; i < size; i++)
    adjacency_rcv[i] = 0;

  /* we build arrays of ghost cell indices for each PE and each level */
  foreach_cell_post (true)
    if (!is_active (cell) && cell.pid != rank) {
      if (cell.neighbors == 0) {
 	free (allocated(0,0));
	allocated(0,0) = NULL;
	point.back->dirty = true;
      }
      else {
	int i;
	for (i = 0; i < nrcv; i++)
	  if (cell.pid == rcv[i])
	    break;
	if (i == nrcv) {
	  rcv = realloc (rcv, ++nrcv*sizeof (int));
	  rcv[nrcv-1] = cell.pid;
	  rcv_halo = realloc (rcv_halo, nrcv*sizeof (CacheLevel *));
	  rcv_halo[nrcv-1] = malloc ((depth() + 1)*sizeof (CacheLevel));
	  for (int j = 0; j <= depth(); j++)
	    cache_level_init (&rcv_halo[nrcv-1][j]);
	}
	cache_level_append (&rcv_halo[i][level], point);
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
  for (int i = 0; i < nrcv; i++)
    for (int j = 0; j <= depth(); j++) {
      CacheLevel halo = rcv_halo[i][j];
      MPI_Request r;
      MPI_Isend (&halo.n, 1, MPI_INT, rcv[i], j, MPI_COMM_WORLD, &r);
      if (halo.n > 0)
	MPI_Isend (halo.p, 2*halo.n, MPI_INT, rcv[i], j, MPI_COMM_WORLD, &r);
    }

  char name[80];
  sprintf (name, "log-%d", rank);
  FILE * logfile = fopen (name, "w");

  /* we receive the ghost cell indices from their respective PEs */
  for (int i = 0; i < size; i++)
    if (adjacency_snd[i]) {
      snd = realloc (snd, ++nsnd*sizeof (int));
      snd[nsnd-1] = i;
      snd_halo = realloc (snd_halo, nsnd*sizeof (CacheLevel *));
      snd_halo[nsnd-1] = malloc ((depth() + 1)*sizeof (CacheLevel));
      for (int j = 0; j <= depth(); j++) {
	CacheLevel * halo = &snd_halo[nsnd-1][j];
	MPI_Status s;
	MPI_Recv (&halo->n, 1, MPI_INT, i, j, MPI_COMM_WORLD, &s);
	if (halo->n > 0) {
	  halo->p = malloc (halo->n*sizeof(IndexLevel));
	  MPI_Recv (halo->p, 2*halo->n, MPI_INT, i, j, MPI_COMM_WORLD, &s);
	  /* we make sure that halo cells are treated as such */
	  foreach_cache(snd_halo[nsnd-1][j], j)
	    cell.neighbors++;
	}
      }
    }

  sprintf (name, "cells-%d", rank);
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  sprintf (name, "ghost-%d", rank);
  fp = fopen (name, "w");
  for (int i = 0; i < nsnd; i++)
    for (int j = 0; j <= depth(); j++)
      foreach_cache (snd_halo[i][j], j)
	fprintf (fp, "%g %g %d %d\n", x, y, snd[i], cell.neighbors);
  fclose (fp);

  scalar s[];
  foreach()
    s[] = 1;

  halo_restriction_mpi (-1, s);
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_cell (b, true) {
      if (is_active (cell) || cell.neighbors > 0) // check this
	s[ghost] = s.boundary[b] (point, s);
      else
	continue;
    }

  foreach_halo_coarse_to_fine (depth())
    if (cell.pid == rank) {
      fprintf (logfile, "%g %g %d %d\n", x, y, cell.pid, cell.neighbors);
      fflush (logfile);
      s[] = (9.*coarse(s,0,0) + 
	     3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
	     coarse(s,child.x,child.y))/16.;  
    }
  for (int l = 0; l <= depth(); l++)
    boundary_mpi (s, l);

#if 0
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_cell (b, true)
      if (!is_active (cell)) { // fixme
	if (cell.neighbors > 0)
	  s[ghost] = s.boundary[b] (point, s);
	else
	  continue;
      }
#endif

  sprintf (name, "a-%d", rank);
  fp = fopen (name, "w");
#if 1
  foreach() {
#if 1
    double lap = (s[1,0] + s[0,1] + s[-1,0] + s[0,-1] - 4.*s[]);
    fprintf (fp, "%g %g %g\n", x, y, lap);
    assert (fabs (lap) < 1e-10);
#else
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	fprintf (fp, "%g %g %g\n", x + i*Delta, y + j*Delta, s[i,j]);
#endif
  }
#else
  foreach_cell()
    if (is_active(cell) || cell.neighbors > 0)
      fprintf (fp, "%g %g %g\n", x, y, s[]);
#endif
  fclose (fp);

  MPI_Finalize();
}
