/**
# Computation of volume fractions from a levelset function */

#include "grid/cartesian.h"
#include "fractions.h"

int main()
{
  X0 = Y0 = -0.5;
  init_grid (32);

  /**
  The interface is a circle centered on the origin and of radius
  0.25. We chose this radius because it leads to degenerate cases where
  the interface intersects the grid exactly on vertices. 
  
  We initialise a levelset function on the vertices of the grid. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(0.25) - sq(x) - sq(y);

  /**
  We then use this function to compute the corresponding volume and
  surface fractions. */

  scalar c[];
  face vector s[];
  fractions (phi, c, s);

  /**
  To check that this is correct, we draw the corresponding facets,
  reconstructed using either `c` and `s` (on stdout), or only `c` (on
  stderr). */

  output_facets (c, stdout, s);
  output_facets (c, stderr);

  free_grid ();
}

/**
This gives this figure where "exact" uses *c* and *s* and "VOF" uses
only *c*.

~~~gnuplot Exact and VOF-reconstucted interface
set output 'plot.png'
set size ratio -1
set key out
plot 'out' w l t "exact", 'log' w l t "VOF"
~~~
*/