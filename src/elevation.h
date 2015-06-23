/**
# Conservation of water surface elevation 

When using the default adaptive reconstruction of variables, the
[Saint-Venant solver](saint-venant.h) will conserve the water depth
when cells are refined or coarsened. However, this will not
necessarily ensure that the "lake-at-rest" condition (i.e. a constant
water surface elevation) is also preserved. In what follows, we
redefine the *refine()* and *coarsen()* methods of the water depth $h$
so that the water surface elevation $\eta$ is conserved.

We start with the reconstruction of fine "wet" cells: */

#if QUADTREE
static void refine_elevation (Point point, scalar h)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  if (h[] >= dry) {
    double eta = zb[] + h[];   // water surface elevation  
    coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
	g.x = gradient (zb[-1,0] + h[-1,0], eta, zb[1,0] + h[1,0])/4.;
    else
      foreach_dimension()
	g.x = (zb[1,0] - zb[-1,0])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta;
      foreach_dimension()
	etac += g.x*child.x;
      h[] = max(0, etac - zb[]);
    }
  }
  else {

    /**
    The "dry" case is a bit more complicated. We look in a 3x3
    neighborhood of the coarse parent cell and compute a depth-weighted
    average of the "wet" surface elevation $\eta$. We need to do this
    because we cannot assume a priori that the surrounding wet cells are
    necessarily close to e.g. $\eta = 0$. */

    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h[i,j] >= dry) {
	  eta += h[i,j]*(zb[i,j] + h[i,j]);
	  v += h[i,j];
	}
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else

      /**
      If none of the surrounding cells is wet, we assume a default sealevel
      at zero. */

      eta = 0.;

    /**
    We then reconstruct the water depth in each child using $\eta$ (of the
    parent cell i.e. a first-order interpolation in contrast to the wet
    case above) and $z_b$ of the child cells. */
    
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta - zb[]);
  }
}

/**
Cell coarsening is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */

static void coarsen_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }

  /**
  ... and use this in combination with $z_b$ (of the coarse cell) to
  compute the water depth $h$.  */
    
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h.refine  = refine_elevation;
  h.coarsen = coarsen_elevation;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif

/**
# "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (*ref*). */

#define radiation(ref) (sqrt (G*max(h[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

/**
# Tide gauges

An array of *Gauge* structures passed to *output_gauges()* will create
a file (called *name*) for each gauge. Each time *output_gauges()* is
called a line will be appended to the file. The line contains the time
and the value of each scalar in *list* in the (wet) cell containing
*(x,y)*. The *desc* field can be filled with a longer description of
the gauge. */

typedef struct {
  char * name;
  double x, y;
  char * desc;
  FILE * fp;
} Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  for (Gauge * g = gauges; g->name; g++) {
    if (!g->fp) {
      g->fp = fopen (g->name, "w");
      if (g->desc)
	fprintf (g->fp, "%s\n", g->desc);
    }
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    Point point = locate (xp, yp);
    if (point.level >= 0 && h[] > dry) {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " %g", s[]);
      fputc ('\n', g->fp);
      fflush (g->fp);
    }
  }
}
