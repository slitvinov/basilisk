/**
# Speed of elementary operations on different grids

This code is compiled using either the default quadtree implementation
or the 'multigrid' regular Cartesian grid implementation.

We use a square, regular, Cartesian grid and vary the resolution from
16^2^ to 2048^2^ (to check for the influence of memory caching). */

#include "poisson.h"
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
  for (int j = 0; j < npe(); j++)
    printf (" %g", mpi[j]/i);
  printf ("\n");
  MPI_Barrier (MPI_COMM_WORLD);
  trace_event (name);
}

int main (int argc, char * argv[])
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : 8;
  int minlevel = argc > 2 ? atoi(argv[2]) : 1;
  timer t;

  init_grid (1 << minlevel);
  mpi_partitioning();

  foreach()
    a[] = b[] = 0.;
  boundary ({a, b});
  poisson (a, b); // to force allocation of extra fields
  
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
    foreach()
      b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[])/sq(Delta);
    boundary ({b});
  }
  mpi_print (t, nloops, tnc*nloops, "laplacian");
  
  /**
     Something simpler: the sum of `a` over the entire mesh. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = npe();
  double sum = 0.;
  while (i--) {
    sum = 0.;
    foreach(reduction(+:sum))
      sum += b[];
  }
  mpi_print (t, nloops, tnc*nloops, "sum");
  fprintf (stderr, "sum: %g\n", sum);

  /**
     The restriction operator. */
  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = 1; // npe();
  while (i--)
    restriction ({b});
  mpi_print (t, nloops, tnc*nloops, "restriction");

  /**
     And finally the Poisson solver. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = 1;
  TOLERANCE = HUGE;
  while (i--)
    poisson (a, b);
  mpi_print (t, nloops, tnc*nloops, "poisson");

  scalar e[];
  foreach()
    e[] = a[] - cos(pi*x)*cos(pi*y);
  fprintf (stderr, "error: %g\n", normf(e).max);
  assert (normf(e).max < 1e-10);
  //  output_ppm (e, file = "error.png");
  
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
