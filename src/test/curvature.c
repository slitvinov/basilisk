/**
# Curvature of a circular interface

This test evaluates the accuracy of the generalised height-function
curvature calculation. It is similar to the case presented in
[Popinet, 2009](references.bib#popinet2009) (Figure 5). The curvatures
of circles with a randomised position and varying radii are computed
and statistics on the error are gathered and displayed on the graph
below. */

#include "fractions.h"
#include "curvature.h"

/**
The function below is called with different arguments for "coarse" and
"fine" interface resolution in order to minimize the computation
times. The number of randomised positions is *nr*, the radius of the
circle is *R*, the maximum refinement level to consider is *levelmax*,
the number of grid points used to initialise (accurately) the volume
fraction field is *ninit^2^* and the statistics for each level are
stored in the array *n*, while the statistics on which method is used
are stored in *sc*. */

void sample_circles (int nr, double R, int levelmax, int ninit,
		     norm * n, cstats * sc)
{
  while (nr--) {

    /**
    It is important to use an accurate enough initialisation of the
    volume fraction field for the circular interface. To do this we
    use the standard (low order) *fractions()* function on a fine
    enough grid. */

    init_grid (ninit);
    scalar c[], kappa[];
    vertex scalar phi[];
    double x0 = noise()/8., y0 = noise()/8.;
    foreach_vertex()
      phi[] = - (sq(x - x0) + sq(y - y0) - sq(R));
    fractions (phi, c);

    /**
    We then successively coarsen this fine initial grid to compute the
    curvature on coarser and coarser grids (thus saving on the
    expensive initial condition). */

    for (int l = levelmax; l >= 3; l--) {
      unrefine (level >= l, {c});
      cstats s = curvature (c, kappa);

      /**
      We store statistics on the methods used for curvature
      computation... */

      sc[l].h += s.h; sc[l].f += s.f; sc[l].a += s.a; sc[l].c += s.c;
      foreach()
	if (c[] > 0. && c[] < 1.) {

          /**
	  ...and error statistics (for a given level of refinement *l*). */
           
	  double e = fabs(kappa[] - 1./R)*R;
	  n[l].volume += dv();
	  n[l].avg += dv()*e;
	  n[l].rms += dv()*e*e;
	  if (e > n[l].max)
	    n[l].max = e;
	}
    }
  }  
}

int main()
{
  origin (-0.5, -0.5);
 
  /**
  We try a wide enough range of radii. */

  for (double R = 0.1; R <= 0.2; R *= 1.2) {

    /**
    We initialize the arrays required to store the statistics for each
    level of refinement. */

    int levelmax = 7;
    norm n[levelmax + 1];
    cstats sc[levelmax + 1];
    for (int i = 0; i <= levelmax; i++) {
      n[i].volume = n[i].avg = n[i].rms = n[i].max = 0;
      sc[i].h = sc[i].f = sc[i].a = sc[i].c = 0.;
    }

    /**
    The finer the (initialisation) mesh, the more expensive the
    computation. On the other hand, we can limit randomisation for the
    higher resolutions (since we expect less "special cases" on fine
    meshes). We thus limit the total runtime by sampling many (100)
    locations on coarse meshes but only few (1) location on the finest
    mesh. Note that to get second-order convergence, the finest mesh
    (128^2^) requires an initialisation at 2048^2^. */

    sample_circles (100, R, 4, 256, n, sc);
    sample_circles (10, R, 6, 512, n, sc);
    sample_circles (1, R, levelmax, 2048, n, sc);

    /**
    Finally we output the statistics for this particular radius and for
    each level of refinement. */

    for (int l = levelmax; l >= 3; l--) {
      n[l].avg /= n[l].volume;
      n[l].rms = sqrt(n[l].rms/n[l].volume);
      double t = sc[l].h + sc[l].f + sc[l].a + sc[l].c;
      fprintf (stderr, "%g %g %g %g %g %g %g %g\n",
	       2.*R*(1 << l),
	       n[l].avg, n[l].rms, n[l].max,
	       sc[l].h/t, sc[l].f/t, sc[l].a/t, sc[l].c/t);
    }
  }

  /**
  At the end of the run, we sort the data by increasing order of
  diameter (in grid points). */

  fflush (stderr);
  system ("sort -k1,2 -n log > log.s; mv -f log.s log");
}

/**
The results are summarised in the figure below. There are two sets of
points: error norms (bottom three curves) and percentages of
curvatures computed with each method (top four curves). As expected
the second-order convergence of the max and RMS norms is recovered for
the pure HF method, for a diameter greater than about 15 grid
points. The RMS norm shows consistent second-order convergence across
the whole range of diameters. There is a significant degradation of
the max norm between 7 and 15 grid points, corresponding with the
introduction of the "averaging" and "(mixed) HF fit" methods. This
degradation is significantly less pronounced for the method
implemented in [Popinet, 2009](references.bib#popinet2009).

The grading of the methods used for curvature calculation follows what
is expected: 

* exclusively HF for $D > 15$, 
* a combination of HF, nearest-neighbor average and "mixed HF fit" for
$3 < D < 15$
* and exclusively "centroids fit" for $D < 3$.

~~~gnuplot Relative curvature error as a function of resolution
set logscale
set grid
set key bottom left
set xlabel 'Diameter (grid points)'
set ylabel 'Relative curvature error / percentage' 
f(x)=(x > 0. ? 100.*x : 1e1000)
plot 2./(x*x) t '2/x^{2}', 'log' u 1:4 w lp t 'Max', '' u 1:3 w lp t 'RMS', \
  '../popinet.csv' u ($1*2):2 w lp t 'Popinet (2009)', \
  'log' u 1:(f($5)) w lp t 'HF', '' u 1:(f($6)) w lp t 'HF fit', \
  '' u 1:(f($7)) w lp t 'Average', '' u 1:(f($8)) w lp t 'Centroids'
~~~
*/
