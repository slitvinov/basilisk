/**
# Axisymmetric coordinates

For problems with a symmetry of revolution around the $z$-axis of a
[cylindrical coordinate
system](http://en.wikipedia.org/wiki/Cylindrical_coordinate_system). The
longitudinal coordinate ($z$-axis) is *x* and the radial coordinate
($\rho$- or $r$-axis) is *y*. Note that *y* (and so *Y0*) cannot be
negative.

We first define a macro which will be used in some geometry-specific
code (e.g. [curvature computation](curvature.h)). */

#define AXI 1

event defaults (i = 0) {

  /**
  By default *cm* is a constant scalar field. To make it variable, we
  need to allocate a new field. */

  if (is_constant(cm))
    cm = new scalar;

  /**
  The volume/area of a cell is proportional to $r$ (i.e. $y$). We need
  to set boundary conditions at the top and bottom so that *cm* is
  interpolated properly when refining/coarsening the mesh. */

  foreach()
    cm[] = y;
  cm[top] = dirichlet(y);
  cm[bottom] = dirichlet(y);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the axis of revolution is zero ($y=r=0$ on the axis). To avoid
  division by zero we set it to epsilon (note that mathematically the
  limit is well posed). */

  if (is_constant(fm.x))
    fm = new face vector;
  foreach_face()
    fm.x[] = max(y, 1e-20);
  fm.t[top] = dirichlet(y);
  fm.t[bottom] = dirichlet(y);

  boundary ({cm, fm});
}
