#define MULTIGRID 1

#include "cartesian-common.h"

@ifndef foreach_level_or_leaf
@ define foreach_level_or_leaf     foreach_level
@ define end_foreach_level_or_leaf end_foreach_level
@endif

@ifndef foreach_boundary_fine_to_coarse
@ def foreach_boundary_fine_to_coarse(dir) {
  for (int _l = depth() - 1; _l >= 0; _l--)
    foreach_boundary_level (dir,_l,false) {
      point.i += ig; point.j += jg;
      ig = -ig; jg = -jg;
      POINT_VARIABLES;
@
@ def end_foreach_boundary_fine_to_coarse()
      ig = -ig; jg = -jg;
    } end_foreach_boundary_level() 
  }
@
@endif

// scalar attributes

attribute {
  double (* prolongation) (Point, scalar);
};

// Multigrid methods

void (* boundary_restriction) (scalar *);

void coarsen_average (Point point, scalar s)
{
  s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
}

void restriction (scalar * list)
{
  scalar * listc = NULL;
  vector * lists = NULL;
  for (scalar s in list) 
    if (!is_constant(s)) {
      if (s.face)
	lists = vectors_add (lists, s.v);
      else
	listc = list_add (listc, s);
    }
  if (lists)
    boundary_normal (lists);
  if (lists || listc)
    foreach_fine_to_coarse() {
      for (scalar s in listc)
	s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
      for (vector v in lists)
	foreach_dimension()
	  v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
    }
  free (listc);
  if (lists) {
    foreach_boundary_fine_to_coarse(right)
      for (vector v in lists)
	v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
    foreach_boundary_fine_to_coarse(top)
      for (vector v in lists)
	v.y[] = (fine(v.y,0,0) + fine(v.y,1,0))/2.;
    free (lists);
  }
}

void wavelet (scalar s, scalar w)
{
  restriction ({s});
  boundary_restriction ({s});
  foreach_fine_to_coarse() {
    if (s.prolongation)
      /* difference between fine value and its prolongation */
      foreach_child()
	w[] = s[] - s.prolongation (point, s);
    else
      /* difference between fine value and bilinearly-interpolated
	 coarse value */
      foreach_child()
	w[] = s[] - (9.*coarse(s,0,0) + 
		     3.*(coarse(s,child.x,0) + coarse(s,0,child.y)) + 
		     coarse(s,child.x,child.y))/16.;
  }
  /* root cell */
  foreach_level(0) w[] = 0.;
}

void refine_bilinear (Point point, scalar s)
{
  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      /* bilinear interpolation from coarser level */
      fine(s,k,l) = (9.*s[] + 
		     3.*(s[2*k-1,0] + s[0,2*l-1]) + 
		     s[2*k-1,2*l-1])/16.;
}

void refine_linear (Point point, scalar s)
{
  struct { double x, y; } g;
  if (s.gradient)
    foreach_dimension()
      g.x = s.gradient (s[-1,0], s[], s[1,0]);
  else
    foreach_dimension()
      g.x = (s[1,0] - s[-1,0])/2.;

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      /* linear interpolation from coarser level (conservative) */
      fine(s,k,l) = s[] + (g.x*(2*k-1) + g.y*(2*l-1))/4.;
}

void refine_reset (Point point, scalar v)
{
  /* foreach_child() */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      fine(v,k,l) = 0.;
}

static void nothing() {}
void * none = nothing;

void coarsen_face (Point point, scalar s)
{
  vector v = s.v;
  foreach_dimension()
    v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
}

void multigrid_boundary_level (scalar * list, int l)
{
  if (l < 0)
    l = depth();
  // traverse the boundaries of a given level
  for (int b = 0; b < nboundary; b++)
    foreach_boundary_level (b, l, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = s.boundary[b] (point, s);
}

void multigrid_boundary_restriction (scalar * list)
{
  // traverse the boundaries of all coarse levels (i.e. not the leaves)
  for (int l = 0; l < depth(); l++)
    multigrid_boundary_level (list, l);
}

void multigrid_debug (Point point)
{
  cartesian_debug (point);

  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
    double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
    for (int k = 0; k <= 1; k++)
      for (int l = 0; l <= 1; l++) {
	for (scalar v in all)
	  fprintf (fp, "%g %g %g ", 
		   xc + k*child.x*Delta*2. + v.d.x*Delta, 
		   yc + l*child.y*Delta*2. + v.d.y*Delta,
		   coarse(v,k*child.x,l*child.y));
	fputc ('\n', fp);
      }
    fclose (fp);
    fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
    double xf = x - Delta/4., yf = y - Delta/4.;
    for (int k = 0; k <= 2; k++)
      for (int l = 0; l <= 2; l++) {
	for (scalar v in all)
	  fprintf (fp, "%g %g %g ", 
		   xf + k*Delta/2. + v.d.x*Delta/4., 
		   yf + l*Delta/2. + v.d.y*Delta/4.,
		   fine(v,k,l));
	fputc ('\n', fp);
      }
    fclose (fp);
    fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  boundary_level = multigrid_boundary_level;
  boundary_restriction = multigrid_boundary_restriction;
}
