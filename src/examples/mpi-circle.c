/**
# Parallel efficiency for a Poisson problem

We run the program below in parallel using

~~~bash
for i in `seq 1 1 8`; do
    mpirun -np $i ./mpi-circle 2> log-$i > out-$i
    sed "s/# /$i pe, /g" < out-$i
done > out-all
mpicc -std=c99 -D_GNU_SOURCE=1 -D_MPI=1 -DFAKE_MPI=1 -O2 -Wall _mpi-circle.c -o mpi-circle -lm
... > out-fake
~~~

The mesh and solution is the same as for [/src/test/circle.c]().

The number of leaf points per core is $N/p$ (assuming perfect
balancing), with $p$ the number of cores and $N$ the total number of
leaves.

The number of points to transfer with MPI will scale like the length
of the boundary i.e. $p^{-1/2}$.

This is confirmed by the data

~~~gnuplot Number of points as a function of the number of cores
set logscale
set xlabel 'Number of cores'
set ylabel 'Number of points'
N=114880
f(x)=a*x**-0.5
fit [4:] f(x) '< grep MPI out-all' u 1:11 via a
set xtics 1,2,8
plot '< grep MPI out-all' u 1:5 w lp t 'leaves', N/x t 'N/p', \
     '< grep MPI out-all' u 1:11 w lp t 'MPI', f(x) t 'p^{-1/2}'
~~~

A model for the wall-clock time $T$ is
$$
T_p = T_1 p^{-1} + M p^{-1/2}
$$
where $T_1$ is the time taken on a single core and $M$ is a measure of
the performance of MPI transfers. This gives

~~~gnuplot wall-clock time as a function of the number of cores
set ylabel 'Wall-clock time'
T(p)=T1/p+M/sqrt(p)
fit T(x) '< grep Quadtree out-all' u 1:8 via T1,M
f(x)=a/x
fit f(x) '< grep Quadtree out-all' u 1:8 via a
f1(x)=b/sqrt(x)
fit f1(x) '< grep Quadtree out-all' u 1:8 via b
set yrange [:20]
plot '< grep Quadtree out-all' u 1:8 w lp t '', T(x), \
     f(x) t 'p^{-1}', f1(x) t 'p^{-1/2}'
~~~

If $M$ is not much smaller than $T_1$, the time will scale like
$p^{-1/2}$ rather than the optimal $p^-1$. This model assumes that the
"serial speed" $T_1$ is independent from $p$. In practice this is not
the case as serial jobs will compete for cache etc...

To test this, we can plot the computing time minus the time spent in
MPI calls, this should scale like $T_1/p$, while the time spent in MPI
calls should scale like $M p^{-1/2}$.

~~~gnuplot Compute time and MPI time
fit [0:2] f(x) "< awk '/Quadtree/{pe=$1;t=$8}/MPI/{print pe,t-$7}' < out-all" via a
set yrange [*:*]
unset logscale
plot									\
  '< grep Quadtree out-all' u 1:8 w lp t 'Total',			\
  '< grep MPI out-all' u 1:7 w lp t 'MPI',				\
  "< awk '/Quadtree/{pe=$1;t=$8}/MPI/{print pe,t-$7}' < out-all" w lp	\
  t 'Total - MPI', f(x) t 'p^{-1}'
~~~

This graph shows that the overhead comes mostly from the serial part
of the code which slows down as the number of cores increases. 

To check whether this is due to interactions between the serial part
of the code and MPI (e.g. the serial code and MPI competing for cache
etc...), we rerun but with FAKE_MPI set, this removes all MPI calls
(but still uses the same partitioning). We get

~~~gnuplot Compute time and MPI time (FAKE_MPI)
fit [0:2] f(x) "< awk '/Quadtree/{pe=$1;t=$8}/MPI/{print pe,t-$7}' < out-fake" via a
plot									\
  '< grep Quadtree out-fake' u 1:8 w lp t 'Total',			\
  '< grep MPI out-fake' u 1:7 w lp t 'MPI',				\
  "< awk '/Quadtree/{pe=$1;t=$8}/MPI/{print pe,t-$7}' < out-fake" w lp	\
  t 'Total - MPI', f(x) t 'p^{-1}'
~~~

This shows that the effect of MPI on the serial part is negligible.

We can also compare the MPI implementation with the OpenMP version, using

~~~bash
gcc -std=c99 -D_GNU_SOURCE=1 -fopenmp -O2 -Wall _mpi-circle.c -o mpi-circle -lm
for i in `seq 1 1 8`; do
    OMP_NUM_THREADS=$i ./mpi-circle 2> log-$i > out-$i
    sed "s/# /$i pe, /g" < out-$i
done > out-omp
~~~

We get the comparison

~~~gnuplot Comparison between MPI and OpenMP
fit [0:2] f(x) 'out-omp' u 1:8 via a
set yrange [:20]
set logscale
plot									\
  '< grep Quadtree out-all' u 1:8 w lp t 'MPI',				\
  'out-omp' u 1:8 w lp t 'OpenMP',					\
  f(x) t 'p^{-1}'
~~~
*/

#include "utils.h"
#include "poisson.h"

double solution (double x, double y)
{
  return sin(3.*pi*x)*sin(3.*pi*y);
}

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

void solve (int depth)
{
  origin (-0.5, -0.5);
  init_grid(1);

  while (refine_function (refine_circle, &depth, NULL));
@if _MPI
  mpi_partitioning();
@endif

  scalar a[], b[];

  /* Dirichlet condition on all boundaries */
  a[right]  = dirichlet (solution(x, y));
  a[left]   = dirichlet (solution(x, y));
  a[top]    = dirichlet (solution(x, y));
  a[bottom] = dirichlet (solution(x, y));

  timer start = timer_start();

  foreach() {
    a[] = 0.;
    b[] = -18.*pi*pi*sin(3.*pi*x)*sin(3.*pi*y);
  }
  boundary ({a});

  TOLERANCE = 0.;
  NITERMAX = 100;
  mgstats s = poisson (a, b);

  double max = 0;
  foreach(reduction(max:max)) {
    double e = a[] - solution(x, y);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  
  timer_print (start, s.i, -1);
  fprintf (stderr, "max error %d %g\n", depth, max);

@if _MPI
  mpi_boundary_stats (stdout);
@endif

  free_grid();
}

int main()
{
  solve (10);
}
