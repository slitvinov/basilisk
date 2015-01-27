/**
# Speed of elementary operations on different grids

This code is compiled using either the default quadtree implementation
or the 'multigrid' regular Cartesian grid implementation.

We use a square, regular, Cartesian grid and vary the resolution from
16^2^ to 2048^2^ (to check for the influence of memory caching). */

#include "utils.h"

scalar a[], b[];

static void mpi_print (timer t, int i, size_t tnc,
		       const char * name)
{
  trace_event (name);
  double mpi[npe()];
  timing s = timer_timing (t, i, tnc, mpi);
  printf ("%d %g %g %g %s %g %g %g %.2g%% %ld %ld ",
	  npe(), s.cpu/i, s.real/i, s.speed, name, s.min/i, s.avg/i, s.max/i,
	  100.*s.avg/s.real, s.mem, s.tnc);
#if 0
  for (int j = 0; j < npe(); j++)
    printf (" %g", mpi[j]/i);
#endif
  printf ("\n");
  MPI_Barrier (MPI_COMM_WORLD);
  trace_event (name);
}

static void laplacian()
{
  foreach()
    b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[])/sq(Delta);
  boundary ({b});
}

static void laplacian_overlapping (int maxlevel)
{
  foreach_cache (quadtree->bleaves,)
    b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[])/sq(Delta);
  snd_rcv_send (&((MpiBoundary *)mpi_boundary)->prolongation, {b}, maxlevel);
  foreach_cache (quadtree->ileaves,)
    b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[])/sq(Delta);
  snd_rcv_receive (&((MpiBoundary *)mpi_boundary)->prolongation, {b}, maxlevel);
}

int main (int argc, char * argv[])
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : 8;
  int minlevel = argc > 2 ? atoi(argv[2]) : 1;
  bool overlapping = argc > 3 ? atoi(argv[3]) : false;
  timer t;

  init_grid (1 << minlevel);
  mpi_partitioning();

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  refine (level < maxlevel, NULL);
  size_t tnc = 0;
  foreach(reduction(+:tnc))
    tnc++;
  mpi_print (t, 1, tnc, "refine");

  int nloops, i;

  /**
     We fill `a` with a simple function. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = npe();
  while (i--) {
    foreach()
      a[] = cos(pi*x)*cos(pi*y);
#if 0
    boundary ({a});
#else
    prof_start ("barrier");
    MPI_Barrier (MPI_COMM_WORLD);
    prof_stop();
#endif
  }
  mpi_print (t, nloops, tnc*nloops, "cos");
  boundary ({a});

  /**
     We set a number of loops proportional to the number of grid
     points so that times are comparable on all grids. */

  /**
     Here we compute
     $$
     b = \nabla^2 a
     $$
     using a 5-points Laplacian operator. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = npe();
  while (i--) {
    if (overlapping)
      laplacian_overlapping (maxlevel);
    else
      laplacian();
  }
  mpi_print (t, nloops, tnc*nloops, "laplacian");
    
  int n = 0;
  foreach()
    n++;
  int nmin = n, nmax = n;
  mpi_all_reduce (nmin, MPI_INT, MPI_MIN);
  mpi_all_reduce (nmax, MPI_INT, MPI_MAX);
  fprintf (stderr, "balance %d %d\n", nmin, nmax);
}

/**
## Results

This graph shows the speed of the quadtree implementation for each
operation relative to the speed on the Cartesian mesh. As expected,
the overhead is relatively larger for the simpler operations (e.g. sum
of all elements). This graph is quite sensitive to the exact machine
architecture (cache hierarchy etc...).

![Relative speed of simple operations on a quadtree mesh](laplacian/plot.png)

The absolute speed for the Laplacian operator on both grid
implementations is shown below. Note that Cartesian meshes are fast!
(i.e. hundreds of million of grid points per second).

![Absolute speed of the 5-points Laplacian on Cartesian and quadtree
 meshes](laplacian/speed.png) */
