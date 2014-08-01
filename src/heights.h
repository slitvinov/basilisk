#if !QUADTREE
void heights (scalar c, vector h)
{
  foreach()
    foreach_dimension()
      h.x[] = Delta*(c[-2,0] + c[-1,0] + c[] + c[1,0] + c[2,0] - 7./2.);
  for (int j = -1; j <= 1; j += 2) {
    vector cs[];
    // shift by 2
    foreach()
      foreach_dimension()
        cs.x[] = c[2*j,0];
    // boundary conditions
    boundary ((scalar *){cs});
    // add 3rd cell
    foreach()
      foreach_dimension() {
        h.x[] += Delta*cs.x[j,0];
	if (j == 1 && cs.x[1,0] > 0.5)
	  h.x[] = -h.x[];
      }
  }
  boundary ((scalar *){h});
}
#else // QUADTREE
foreach_dimension()
static void h_refine_x (Point point, scalar h)
{
  for (int k = 0; k <= 1; k++)
    for (int l = 0; l <= 1; l++) {
      // third-order taylor series
      int m = 2*k-1, n = 2*l-1;
      double hf = (15.*h[]/2. + h[0,1]*(1./4. + n) + h[0,-1]*(1./4. - n))/8.;
      fine(h,k,l) = hf - m*Delta/4.;
    }
}

static void coarsen_curvature (Point point, scalar kappa)
{
  double k = 0., s = 0.;
  for (int i = 0; i <= 1; i++)
    for (int j = 0; j <= 1; j++)
      if (fine(kappa,i,j) != nodata) {
	k += fine(kappa,i,j);
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
        h.x[] = Delta*(c[-2,0] + c[-1,0] + c[] + c[1,0] + c[2,0] - 7./2.);
    for (int j = -1; j <= 1; j += 2) {
      vector cs[];
      // shift fine level by 2
      foreach_level (l)
	foreach_dimension()
	  cs.x[] = c[2*j,0];
      // boundary conditions
      // a) shift coarse level by 1 for pairs of adjacent cells
      foreach_level (l-1) 
	foreach_dimension() {
	  cs.x[] = c[j,0];
	  cs.x[j,0] = c[2*j,0];
        }
      // b) prolongate to fine level
      foreach_halo (prolongation, l-1)
	foreach_dimension()
	  c.prolongation (point, cs.x);
      boundary_iterate (halo_prolongation, (scalar *){cs}, l, l);
      // add 3rd cell
      foreach_level (l)
        foreach_dimension() {
	  h.x[] += Delta*cs.x[j,0];
	  if (j == 1 && cs.x[1,0] > 0.5)
	    h.x[] = -h.x[];
        }
      }
  }
  // boundary conditions for h
  foreach_dimension() {
    h.x.prolongation = h_refine_x;
    h.x.coarsen = none;
  }
  boundary ((scalar *){h});
}
#endif // QUADTREE

foreach_dimension()
static double kappa_x (Point point, vector h) {
  double hx = (h.y[1,0] - h.y[-1,0])/(2.*Delta);
  double hxx = (h.y[1,0] + h.y[-1,0] - 2.*h.y[])/sq(Delta);
  return hxx/pow(1. + sq(hx), 3/2.);
}

void curvature (scalar c, scalar kappa)
{
  vector h[];
  heights (c, h);
  foreach() {
    if (c[] > 0. && c[] < 1.) {
      coord n = {c[1,0] - c[-1,0], c[0,1] - c[0,-1]};
      kappa[] = (fabs(n.x) > fabs(n.y) ? 
		 sign(n.x)*kappa_y (point, h) : 
		 sign(n.y)*kappa_x (point, h));
    }
    else
      kappa[] = nodata;
  }

#if QUADTREE
  kappa.prolongation = refine_injection;
  kappa.coarsen = coarsen_curvature;
#endif
  boundary ({kappa});

  // "diffuse" curvature around the interface
  foreach()
    if (kappa[] == nodata) {
      double k = 0., s = 0.;
      for (int i = -1; i <= 1; i++)
	for (int j = -1; j <= 1; j++)
	  if (c[i,j] > 0. && c[i,j] < 1.) {
	    k += kappa[i,j];
	    s++;
	  }
      kappa[] = s ? k/s : nodata;
    }
  boundary ({kappa});
}
