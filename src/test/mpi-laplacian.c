/**
# Parallel scalability

This code can be used to test the scalability of common operations
(traversal, restriction, boundary conditions etc...) and their
combinations, in particular the multigrid Poisson solver. */

#include "poisson.h"
#include "utils.h"

scalar a[], b[];

static void mpi_print (timer t, int i, size_t tnc,
		       const char * name)
{
  double mpi[npe()];
  timing s = timer_timing (t, i, tnc, mpi);

#if 0
  scalar wt[];
  foreach()
    wt[] = mpi[pid()];
  char fname[80];
  sprintf (fname, "%s-%d.ppm", name, npe());
  FILE * fp = pid() ? NULL : fopen (fname, "w");
  output_ppm (wt, fp, n = 512);
#endif
  
  if (pid() == 0) {
    /**
    *s.min/i*, *s.avg/i*, *s.max/i* are the minimum, average and maximum
    *times spent in communication routines. */

    printf ("%d %g %g %g %s %g %g %g %.2g%% %ld %ld ",
	    npe(), s.cpu/i, s.real/i, s.speed, name, s.min/i, s.avg/i, s.max/i,
	    100.*s.avg/s.real, s.mem, s.tnc);

    /**
    We also output the times spent in communication routines for each process. */
#if 0
    for (int j = 0; j < npe(); j++)
      printf (" %g", mpi[j]/i);
#endif
    printf ("\n");
  }
  
  MPI_Barrier (MPI_COMM_WORLD);
}

int main (int argc, char * argv[])
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : (dimension == 2 ? 8 : 5);
  int minlevel = argc > 2 ? atoi(argv[2]) : 1;
  timer t;

  init_grid (1 << minlevel);

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
      a[] = cos(pi*x)*cos(pi*y)*cos(pi*z);
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
  Here we compute
  $$
  b = \nabla^2 a
  $$
  using a 5-points Laplacian operator. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = npe();
  while (i--) {
    foreach() {
      b[] = 0.;
      foreach_dimension()
        b[] += a[1] + a[-1];
      b[] = (b[] - 2.*dimension*a[])/sq(Delta);
    }
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
  if (pid() == 0)
    fprintf (stderr, "sum: %g\n", sum);

  /**
  The restriction operator. */

  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  i = nloops = 1;
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
    e[] = a[] - cos(pi*x)*cos(pi*y)*cos(pi*z);
  double max = normf(e).max;
  if (pid() == 0)
    fprintf (stderr, "error: %g\n", max);
  assert (normf(e).max < 1e-10);
  //  output_ppm (e, file = "error.png");

  sum = 0.;
  int n = 0;
  foreach() {
    e[] = (long) &(a[]);
    sum += e[];
    n++;
  }
  foreach()
    e[] -= sum/n;
#if 0
  FILE * fp = pid() ? NULL : fopen ("map.ppm", "w");
  output_ppm (e, fp, n = 512);
#endif
  
  int nmin = n, nmax = n;
  mpi_all_reduce (nmin, MPI_INT, MPI_MIN);
  mpi_all_reduce (nmax, MPI_INT, MPI_MAX);
  if (pid() == 0)
    fprintf (stderr, "balance %d %d\n", nmin, nmax);
}

/**
## How to run on Curie

This test is run on
[Curie](http://www.prace-ri.eu/best-practice-guide-curie-html/#curie-configuration). The
C99 source code is generated on a system with *qcc* installed and then
copied to Curie using something like

~~~bash
CC99='cc -std=c99' qcc -source mpi-laplacian.c
scp _mpi-laplacian.c popinets@curie-fr.ccc.cea.fr:
~~~

On Curie the following script (*run.sh*) is used to compile and run
the code

~~~bash
#!/bin/bash
#MSUB -r mpi-laplacian
#MSUB -n 32
#MSUB -T 600
####MSUB -Q test
#MSUB -o basilisk_%I.out
#MSUB -e basilisk_%I.log
#MSUB -q standard
#MSUB -A gen7325
#MSUB -w

LEVEL=10

set -x
cd ${BRIDGE_MSUB_PWD}
mpicc -O2 -Wall -std=c99 -D_MPI=1 _mpi-laplacian.c -o mpi-laplacian -lm
rm -f trace-*
ccc_mprun -n ${BRIDGE_MSUB_NPROC} ./mpi-laplacian ${LEVEL} 8 \
    2> log-${LEVEL}-${BRIDGE_MSUB_NPROC} > out-${LEVEL}-${BRIDGE_MSUB_NPROC}
~~~

LEVEL is varied from 10 (~1 million grid points) to 15 (~1 billion
grid points) and the number of processes up to 16384 using something like

~~~bash
for i in 16 32 64 128 256 512 1024 2048 4096 8192 16384; do 
  ccc_msub -n $i run.sh; 
done
~~~

The results can then be collected using

~~~bash
tar czvf curie.tgz out-*-*
~~~

## Results on Curie

The memory usage per core is given below. The curves are a model (see
[mpi-laplacian.plot]() for details). The increase for a large number
of cores corresponds to the memory overhead of communication buffers
etc...

![Memory usage on Curie](curie/memory.png)

The wall-clock time for one iteration of the multigrid Poisson solver
is given below. The red lines are a model of strong scaling. For a low
enough number of cores, close to perfect scaling is obtained with a
best fit computation speed close to 2.1 million grid points/core.

The pink lines connect points corresponding with weak (or
*iso-granular*) scaling i.e. quadrupling both the computation size and
the number of cores. The ideal weak scaling would give horizontal
lines (i.e. constant computation time for proportionally-increasing
problem sizes and number of cores).

![Wall-clock time on Curie for the Poisson solver](curie/poisson.png)

The time spent in communication routines is illustrated below,
together with the model (corresponding to the communication part of
the total time in the previous figure).

![Communication time on Curie for the Poisson solver](curie/poisson-mpi.png)

Similar results are obtained for a pure Laplacian with a best fit
speed of order 30 million grid points/core. This speed is larger than
for the Poisson solver and hence does not scale as well.

![Wall-clock time on Curie for the Laplacian](curie/laplacian.png)

![Communication time on Curie for the Laplacian](curie/laplacian-mpi.png)

Similarly, the pure restriction operator scales quite well, with a
best fit speed of around 35 million grid points/core.

![Wall-clock time on Curie for the restriction operator](curie/restriction.png)

![Communication time on Curie for the restriction operator](curie/restriction-mpi.png)

*/
