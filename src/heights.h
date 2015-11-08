/**
# Height-Functions

The "height-function" is a vector field which gives the distance,
along each coordinate axis, from the center of the cell to the closest
interface defined by a volume fraction field. This distance is
estimated using the "column integral" of the volume fraction in the
corresponding direction. This integral is not always defined (for
example because the interface is too far i.e. farther than 3.5 cells
in our implementation) in which case the value of the field is set to
*nodata*. See e.g. [Popinet, 2009](references.bib#popinet2009) for
more details on height functions.

We also store the "orientation" of the height function together with its
value by adding 10 if the volume fraction is unity on the "top" end. The
function below applied to the value will return the corresponding height
and orientation. 

The distance is normalised with the cell size so that the coordinates
of the interface are given by

~~~c
(x, y + Delta*height(h.y[])) or (x + Delta*height(h.x[]), y)
~~~
*/

static inline double height (double H) {
  return H > 5. ? H - 10. : H < -5. ? H + 10. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 5.;
}

/**
## Half-column integration 

This helper function performs the integration on half a column, either
"downward" (*j = -1*) or "upward" (*j = 1*). */

static void half_column (Point point, scalar c, vector h, vector cs, int j)
{
  /**
  We store both the state and the temporary value of the height
  function in two floats stored in the (double) height function
  field. The 'state' of the height function can be: *complete* if both
  ends were found, zero or one if one end was found and between zero
  and one if only the interface was found. */

  const int complete = -1;

  foreach_dimension() {

    /**
     *S* is the state and *H* the (partial) value of the height
     function. If we are on the (first) downward integration (*j =
     -1*) we initialise *S* and *H* with the volume fraction in
     the current cell. */

    double S = c[], H = S, ci, a;

    /**
     On the upward integration (*j = 1*), we recover the state of the
     downward integration. Both the state and the (possibly partial)
     height value are encoded in a single number using a base 100
     shift for the state. */
    
    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {
      
      /**
      Check whether this is an inconsistent height. */

      if (h.x[] == 300.)
	state.s = complete, state.h = nodata;

      /**
      Otherwise, this is either a complete or a partial height. */
      
      else {
	int s = (h.x[] + 5.)/100.;
	state.h = h.x[] - 100.*s;
	state.s = s - 1;
      }

      /**
      If this is a complete height, we start a fresh upward
      integration. */
      
      if (state.s != complete)
	S = state.s, H = state.h;
    }
    
    /**
     We consider the four neighboring cells of the half column, the
     corresponding volume fraction *ci* is recovered either from the
     standard volume fraction field *c* (first two cells) or from the
     shifted field *cs* (last two cells). The construction of *cs* is
     explained in the next section. */
    
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? c[i*j] : cs.x[(i - 2)*j];
      H += ci;
      
      /**
       We then check whether the partial height is complete or not. */
      
      if (S > 0. && S < 1.) {
	S = ci;
	if (ci <= 0. || ci >= 1.) {
	  
	  /**
	   We just left an interfacial cell (*S*) and found a full or
	   empty cell (*ci*): this is a partial height and we can stop
	   the integration. If the cell is full (*ci = 1*) we shift
	   the origin of the height. */
	  
	  H -= i*ci;
	  break;
	}
      }
      
      /**
       If *S* is empty or full and *ci* is full or empty, we went
       right through he interface i.e. the height is complete and
       we can stop the integration. The origin is shifted
       appropriately and the orientation is encoded using the "10
       trick". */
      
      else if (S >= 1. && ci <= 0.) {
	H = (H - 0.5)*j + (j == -1)*10.;
	S = complete;
	break;
      }
      else if (S <= 0. && ci >= 1.) {
	H = (i + 0.5 - H)*j + (j == 1)*10.;
	S = complete;
	break;
      }
      
      /**
       If *ci* is identical to *S* (which is empty or full), we
       check that *H* is an integer (i.e. its fractional value is
       zero), otherwise we are in the case where we found an
       interface but didn't go through it: this is an
       inconsistent height and we stop the integration. */
      
      else if (S == ci && modf(H, &a))
	break;
    }

    /**
     We update the global state using the state *S* of the
     half-integration. */
    
    if (j == -1) {

      /**
       For the downward integration, we check that the partial heights
       (*S != complete*) are consistent: if the first cell is full
       or empty or if the last cell is interfacial, the partial
       height is marked as inconsistent. */
      
      if (S != complete && ((c[] <= 0. || c[] >= 1.) ||
			    (S > 0. && S < 1.)))
	h.x[] = 300.; // inconsistent
      else if (S == complete)
	h.x[] = H;
      else

	/**
	This is a partial height: we encode the state using a base 100
	shift. */
	
	h.x[] = H + 100.*(1. + (S >= 1.));
    }
    else { // j = 1

      /**
       For the upward integration, we update the current *state*
       using the state of the half-integration *S* only if the
       first downward integration returned a partial height, or if
       the upward integration returned a complete height with a
       smaller value than the (complete) height of the downward
       integration. */
	  
      if (state.s != complete ||
	  (S == complete && fabs(height(H)) < fabs(height(state.h))))
	state.s = S, state.h = H;
      
      /**
       Finally, we set the vector field *h* using the state and
       height. */
      
      if (state.s != complete)
	h.x[] = nodata;
      else
	h.x[] = (state.h > 1e10 ? nodata : state.h);
    }
  }
}

/**
## Multigrid implementation

The *heights()* function takes a volume fraction field *c* and returns
the height function vector field *h*. */

#if !QUADTREE
trace
void heights (scalar c, vector h)
{

  /**
  We need a 9-points-high stencil (rather than the default
  5-points). To do this we store in *cs* the volume fraction field *c*
  shifted by 2 grid points in the respective directions. We make sure
  that this field uses the same boundary conditions as *c*. */
  
  vector cs[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      cs.x.boundary[i] = c.boundary[i];

  /**
  To compute the height function, we sum the volume fractions in a
  (half-)column starting at the current cell. We start by integrating
  downward (*j = -1*) and then integrate upward (*j = 1*). */
  
  for (int j = -1; j <= 1; j += 2) {

    /**
    We first build the shifted (by $\pm 2$) volume fraction field in each
    direction. */
    
    foreach()
      foreach_dimension()
        cs.x[] = c[2*j];
    boundary ((scalar *){cs});

    /**
    We sum the half-column, downward or upward. */
    
    foreach()
      half_column (point, c, h, cs, j);
  }
  boundary ((scalar *){h});
}

/**
## Quadtree implementation 

We first define the prolongation functions for heights. We first check
that the (three in 2D, nine in 3D) coarse heights are defined and have
compatible orientations. If not, the children heights are
undefined. Otherwise, a parabolic fit of the coarse heights is used to
compute the children heights. */

#else // QUADTREE
foreach_dimension()
static void refine_h_x (Point point, scalar h)
{
#if dimension == 2
  int ori = orientation(h[]);
  for (int i = -1; i <= 1; i++)
    if (h[0,i] == nodata || orientation(h[0,i]) != ori) {
      foreach_child()
        h[] = nodata;
      return;
    }

  double h0 = (30.*height(h[]) + height(h[0,1]) + height(h[0,-1]))/16.;
  double dh = (height(h[0,1]) - height(h[0,-1]))/4.;
  foreach_child()
    h[] = h0 + dh*child.y - child.x/2. + 10.*ori;
#else
  assert (false);
#endif
}

/**
The *heights()* function implementation is similar to the multigrid
case, but the construction of the shifted volume fraction field *cs*
is more complex. */

void heights (scalar c, vector h)
{
  vector cs[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      cs.x.boundary[i] = c.boundary[i];

  /**
  To compute the shifted field, we first need to *restrict* the volume
  fraction on all levels. */
  
  restriction ({c});
  for (int j = -1; j <= 1; j += 2)

    /**
    We traverse the quadtree level by level, from coarse to fine. */
    
    for (int l = 1; l <= depth(); l++) {

      /**
      We construct the ($\pm 2$) shifted field at this level. */
      
      foreach_level (l)
	foreach_dimension()
	  cs.x[] = c[2*j];

      /**
      We then need to apply boundary conditions on the shifted
      field. This is more complex than for a constant resolution grid.

      We first construct the ($\pm 1$) shifted field for the
      immediately coarser level. This is done by copying the volume
      fraction field for pairs of adjacent cells. */
      
      foreach_level (l - 1)
	foreach_dimension() {
	  cs.x[] = c[j];
	  cs.x[j] = c[2*j];
        }

      /**
      We can now use this shifted coarse field (which matches the
      shifted fine field) to apply boundary conditions on coarse/fine
      prolongation halos. */
      
      foreach_halo (prolongation, l - 1)
	foreach_dimension()
	  c.prolongation (point, cs.x);
      boundary_iterate (halo_prolongation, (scalar *){cs}, l, l);

      /**
      We can now sum the half-column at this level, downward or upward
      according to *j*. */

      foreach_level (l)
        half_column (point, c, h, cs, j);
    }

  /**
  We set the prolongation function for *h*. The restriction function
  does nothing as we have already defined *h* on all levels. */
  
  foreach_dimension() {
    h.x.prolongation = refine_h_x;
    h.x.coarsen = no_coarsen;
  }
  boundary ((scalar *){h});
}

/**
# Curvature

The curvature field is defined only in interfacial cells. In all the
other cells it takes the value *nodata*. On quadtrees, we need to
redefine the restriction function to take this into account i.e. the
curvature of the parent cell is the average of the curvatures in the
interfacial child cells. */

static void coarsen_curvature (Point point, scalar kappa)
{
  double k = 0., s = 0.;
  foreach_child()
    if (kappa[] != nodata) {
      k += kappa[];
      s += 1.;
    }
  kappa[] = s ? k/s : nodata;
}
#endif // QUADTREE

/**
To compute the curvature, we estimate the derivatives of the height
functions in a given direction (*x*, *y* or *z*). We first check that
all the heights are defined and that their orientations are the
same. We then compute the curvature as
$$
\kappa = \frac{h_{xx}}{(1 + h_x^2)^{3/2}}
$$
in two dimensions, or
$$
\kappa = \frac{h_{xx}(1 + h_y^2) + h_{yy}(1 + h_x^2) - 2h_{xy}h_xh_y}
{(1 + h_x^2 + h_y^2)^{3/2}}
$$
in three dimensions. */

#if dimension == 2
foreach_dimension()
static double kappa_y (Point point, vector h) {
  int ori = orientation(h.y[]);
  for (int i = -1; i <= 1; i++)
    if (h.y[i] == nodata || orientation(h.y[i]) != ori)
      return nodata;
  double hx = (h.y[1] - h.y[-1])/2.;
  double hxx = (h.y[1] + h.y[-1] - 2.*h.y[])/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}
#else // dimension == 3
foreach_dimension()
static double kappa_z (Point point, vector h) {
  int ori = orientation(h.z[]);
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (h.z[i,j] == nodata || orientation(h.z[i,j]) != ori)
	return nodata;
  double hx = (h.z[1] - h.z[-1])/2.;
  double hy = (h.z[0,1] - h.z[0,-1])/2.;
  double hxx = (h.z[1] + h.z[-1] - 2.*h.z[])/Delta;
  double hyy = (h.z[0,1] + h.z[0,-1] - 2.*h.z[])/Delta;
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#endif

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */

static double height_curvature (Point point, scalar c, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding HF curvature function *kappa*
  (defined above). */

  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;  
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = kappa_x;

  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);
#endif

  /**
  We try each curvature function in turn and return the first
  well-defined value. */
  
  foreach_dimension() {
    double kappa = n.x.kappa (point, h);
    if (kappa != nodata)
      return n.x.n < 0. ? - kappa : kappa;
  }

  /**
  Could not compute curvature using heights. */
  
  return nodata;
}

/**
The function below computes the mean curvature *kappa* of the
interface defined by the volume fraction *c*. */

trace
void curvature (scalar c, scalar kappa)
{
  vector h[];
  heights (c, h);
  foreach() {

    /**
    We first check whether the cell contains an interface or if the
    gradient of *c* on the faces of the cell is non-zero (this is
    important when computing [surface
    tension](tension.h#surface-tension-term)). */
    
    bool interface = (c[] > 0. && c[] < 1.);
    if (!interface && c[] <= 0)
      for (int i = -1; i <= 1 && !interface; i += 2)
	foreach_dimension()
	  if (c[i] >= 1.)
	    interface = true;

    /**
    If we are not close to an interface, we set $\kappa$ to *nodata*. */

    if (!interface)
      kappa[] = nodata;
    
    else {
      kappa[] = height_curvature (point, c, h);
      if (kappa[] != nodata) {

	/**
	We limit the maximum curvature to $1/\Delta$. */
	
	if (fabs(kappa[]) > 1./Delta)
	  kappa[] = sign(kappa[])/Delta;

	/**
	We add the axisymmetric curvature if necessary. */
      
#if AXI
	double nr, r = y, hx;
	coord n;
	foreach_dimension()
	  n.x = c[1] - c[-1];
	if (fabs(n.x) > fabs(n.y)) {
	  hx = (height(h.x[0,1]) - height(h.x[0,-1]))/2.;
	  nr = sign(n.x)*hx;
	}
	else {
	  r += height(h.y[])*Delta;
	  hx = (height(h.y[1,0]) - height(h.y[-1,0]))/2.;
	  nr = - sign(n.y);
	}
	/* limit the minimum radius to half the grid size */
	double kaxi = nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
	kappa[] += kaxi;
#endif
      }
    }
  }

  /**
  On quadtrees we set the prolongation and restriction functions for
  the curvature. */
  
#if QUADTREE
  kappa.prolongation = refine_injection;
  kappa.coarsen = coarsen_curvature;
#endif

  boundary ({kappa});
}
