#define MULTIGRID_MPI 1

typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;

static double mpi_L0 = nodata;
static coord mpi_o = {nodata, nodata, nodata};

int mpi_dims[dimension];

@define BUF void * // workaround for bug in qcc

foreach_dimension()
static BUF snd_x (int i, int dst, int tag, int l, scalar * list,
		  MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  int nl = (1 << l) + 2*GHOSTS;
  size_t size = pow(nl, dimension - 1)*list_len(list)*GHOSTS*sizeof(double);
  double * buf = malloc (size), * b = buf;
  foreach_slice_x (i, i + GHOSTS, l)
    for (scalar s in list)
      *b++ = s[];
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}

foreach_dimension()
static void rcv_x (int i, int src, int tag, int l, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  int nl = (1 << l) + 2*GHOSTS;
  size_t size = pow(nl, dimension - 1)*list_len(list)*GHOSTS*sizeof(double);
  double buf[size], * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
  foreach_slice_x (i, i + GHOSTS, l)
    for (scalar s in list)
      s[] = *b++;
}

trace
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s))
      list1 = list_add (list1, s);
  if (!list1)
    return;

  prof_start ("mpi_boundary_level");
  
  if (level < 0) level = depth();  
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension() {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);  
    MPI_Request reqs[2];
    void * buf[2];
    int nl = (1 << level) + 2*GHOSTS, nr = 0;
    if ((buf[0] = snd_x (nl - 2*GHOSTS, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left,  1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left,  0, level, list1);
    rcv_x (nl - GHOSTS,   right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    free (buf[0]); free (buf[1]);
  }

  free (list1);

  prof_stop();
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  free (m);
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = calloc (1, sizeof (MpiBoundary));
  MPI_Dims_create (npe(), dimension, mpi_dims);
  MPI_Cart_create (MPI_COMM_WORLD, dimension,
		   mpi_dims, Periods, 0, &m->cartcomm);

  // make sure other boundary conditions are not applied
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension() {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  }

  // rescale the x-axis and shift the origin
  int dims[dimension], periods[dimension], coords[dimension];
  MPI_Cart_get (m->cartcomm, dimension, dims, periods, coords);
  
  if (L0 != mpi_L0)
    L0 /= dims[0];
  
  struct { double * x, * y, * z; } o = {&X0, &Y0, &Z0};
  foreach_dimension() {
    if (*o.x != mpi_o.x)
      *o.x += L0*coords[dir.x];
    else if (L0 != mpi_L0)
      *o.x += (L0 - mpi_L0)*coords[dir.x];
    mpi_o.x = *o.x;
  }

  mpi_L0 = L0;

  // setup boundary methods and add to list of boundary conditions
  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}

struct Dimensions {
  int nx, ny, nz;
};
 
void dimensions (struct Dimensions p)
{
  for (int i = 0; i < dimension; i++)
    mpi_dims[i] = (&p.nx)[i];
}
