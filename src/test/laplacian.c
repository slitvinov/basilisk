/**
# Speed of elementary operations on different grids

This code is compiled using either the default quadtree implementation
or the 'multigrid' regular Cartesian grid implementation.

We use a square, regular, Cartesian grid and vary the resolution from
16^2^ to 2048^2^ (to check for the influence of memory caching). */

#include "utils.h"

scalar a[], b[];

int main ()
{
  for (int l = 4; l <= 11; l++) {
    init_grid (1 << l);
    int nloops, i;
    clock_t start, end;

/**
We fill `a` with a simple function. */

    foreach()
      a[] = cos(2.*pi*x)*sin(2.*pi*y);
    boundary ({a});

/** 
We set a number of loops proportional to the number of grid
points so that times are comparable on all grids. */

    nloops = i = (1 << 25) >> 2*l;

/**
Here we compute
$$
b = \nabla^2 a
$$
using a 5-points Laplacian operator. */

    start = clock();
    while (i--)
      foreach()
	b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[]);
    end = clock();
    fprintf (stderr, "%d %g\n", l, 
	     1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)));

/**
Something simpler: the sum of `a` over the entire mesh. */

    nloops = i = (1 << 25) >> 2*l;
    double sum = 0.;
    start = clock();
    while (i--)
      foreach()
	sum += a[];
    end = clock();
    printf ("sum %d %g %g\n", l, 
	    1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)), sum);

/**
And finally the restriction operator. */

    nloops = i = (1 << 25) >> 2*l;
    start = clock();
    while (i--)
      restriction ({b});
    end = clock();
    printf ("res %d %g %g\n", l, 
	    1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)), sum);

    free_grid();
  }
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
