#if !QUADTREE
void heights (scalar c, vector h)
{
  foreach()
    foreach_dimension()
      h.x[] = Delta*(c[-2] + c[-1] + c[] + c[1] + c[2] - 7./2.);
  for (int j = -1; j <= 1; j += 2) {
    vector cs[];
    // shift by 2
    foreach()
      foreach_dimension()
        cs.x[] = c[2*j];
    // boundary conditions
    boundary ((scalar *){cs});
    // add 3rd cell
    foreach()
      foreach_dimension() {
        h.x[] += Delta*cs.x[j];
	if (j == 1 && cs.x[1] > 0.5)
	  h.x[] = -h.x[];
      }
  }
  boundary ((scalar *){h});
}
#else // QUADTREE
foreach_dimension()
static void h_refine_x (Point point, scalar h)
{
#if dimension == 2
  if (h[] == nodata || h[0,1] == nodata || h[0,-1] == nodata)
    foreach_child()
      h[] = nodata;
  else {
    double h0 = (15.*h[]/2. + (h[0,1] + h[0,-1])/4.)/8.;
    double dh = (h[0,1] - h[0,-1])/8.;
    foreach_child()
      h[] = h0 + dh*child.y - child.x*Delta/2.;
  }
#else
  assert (false);
#endif
}

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

void heights (scalar c, vector h)
{
  restriction ({c});
  for (int l = 1; l <= depth(); l++) {
    foreach_level (l)
      foreach_dimension()
        h.x[] = Delta*(c[-2] + c[-1] + c[] + c[1] + c[2] - 7./2.);
    for (int j = -1; j <= 1; j += 2) {
      vector cs[];
      // shift fine level by 2
      foreach_level (l)
	foreach_dimension()
	  cs.x[] = c[2*j];
      // boundary conditions
      // a) shift coarse level by 1 for pairs of adjacent cells
      foreach_level (l-1)
	foreach_dimension() {
	  cs.x[] = c[j];
	  cs.x[j] = c[2*j];
        }
      // b) prolongate to fine level
      foreach_halo (prolongation, l-1)
	foreach_dimension()
	  c.prolongation (point, cs.x);
      boundary_iterate (halo_prolongation, (scalar *){cs}, l, l);
      // add 3rd cell
      foreach_level (l)
        foreach_dimension() {
	  h.x[] += Delta*cs.x[j];
	  if (j == 1 && cs.x[1] > 0.5)
	    h.x[] = -h.x[];
        }
      }
  }
  // boundary conditions for h
  foreach_dimension() {
    h.x.prolongation = h_refine_x;
    h.x.coarsen = no_coarsen;
  }
  boundary ((scalar *){h});
}
#endif // QUADTREE

#if dimension == 2
foreach_dimension()
static double kappa_y (Point point, vector h) {
  for (int i = -1; i <= 1; i++)
    if (h.y[i] == nodata)
      return nodata;
  double hx = (h.y[1] - h.y[-1])/(2.*Delta);
  double hxx = (h.y[1] + h.y[-1] - 2.*h.y[])/sq(Delta);
  return hxx/pow(1. + sq(hx), 3/2.);
}
#else // dimension == 3
foreach_dimension()
static double kappa_z (Point point, vector h) {
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (h.z[i,j] == nodata)
	return nodata;
  double hx = (h.z[1] - h.z[-1])/(2.*Delta);
  double hy = (h.z[0,1] - h.z[0,-1])/(2.*Delta);
  double hxx = (h.z[1] + h.z[-1] - 2.*h.z[])/sq(Delta);
  double hyy = (h.z[0,1] + h.z[0,-1] - 2.*h.z[])/sq(Delta);
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*sq(Delta));
  return (hxx + hyy + hxx*sq(hy) + hyy*sq(hx) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#endif

void curvature (scalar c, scalar kappa)
{
  vector h[];
  heights (c, h);
  foreach() {
    bool interface = (c[] > 0. && c[] < 1.);
    if (!interface && c[] <= 0)
      for (int i = -1; i <= 1 && !interface; i += 2)
	foreach_dimension()
	  if (c[i] >= 1.)
	    interface = true;
    if (interface) {
      double nmax = 0.;
      kappa[] = nodata;
      foreach_dimension() {
	double n = c[1] - c[-1];
	if (fabs(n) > nmax) {
	  kappa[] = kappa_x (point, h);
	  if (kappa[] != nodata) {
	    nmax = fabs(n);
	    if (n < 0.)
	      kappa[] = - kappa[];
	    if (fabs(kappa[]) > 1./Delta)
	      kappa[] = sign(kappa[])/Delta;
	  }
	}
      }
#if AXI
      double nr, r = y, hx;
      coord n;
      foreach_dimension()
	n.x = c[1] - c[-1];
      if (fabs(n.x) > fabs(n.y)) {
	hx = (h.x[0,1] - h.x[0,-1])/(2.*Delta);
	nr = sign(n.x)*hx;
      }
      else {
	r += h.y[];
	hx = (h.y[1,0] - h.y[-1,0])/(2.*Delta);
	nr = - sign(n.y);
      }
      /* limit the minimum radius to half the grid size */
      double kaxi = nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
      kappa[] += kaxi;
#endif
    }
    else
      kappa[] = nodata;
  }

#if QUADTREE
  kappa.prolongation = refine_injection;
  kappa.coarsen = coarsen_curvature;
#endif
  boundary ({kappa});
}
