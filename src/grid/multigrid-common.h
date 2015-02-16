#define MULTIGRID 1

#include "cartesian-common.h"

@ifndef foreach_level_or_leaf
@ define foreach_level_or_leaf     foreach_level
@ define end_foreach_level_or_leaf end_foreach_level
@endif

@ifndef foreach_coarse_level
@ define foreach_coarse_level      foreach_level
@ define end_foreach_coarse_level  end_foreach_level
@endif

// scalar attributes

attribute {
  void (* prolongation) (Point, scalar);
  void (* coarsen)      (Point, scalar);
}

// Multigrid methods

void coarsen_average (Point point, scalar s)
{
  s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
}

void coarsen_volume_average (Point point, scalar s)
{
  s[] = 0.;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      s[] += fine(cm,i,j)*fine(s,i,j);
  s[] /= 4.*cm[];
}

void restriction (scalar * list)
{
  scalar * listc = NULL;
  vector * listf = NULL;
  for (scalar s in list) 
    if (!is_constant(s)) {
      if (s.face)
	listf = vectors_add (listf, s.v);
      else
	listc = list_add (listc, s);
    }
  if (listf)
    boundary_flux (listf);
  if (listf || listc) {
    boundary_iterate (restriction, list, depth());
    for (int l = depth() - 1; l >= 0; l--) {
      foreach_coarse_level(l) {
	// fixme: this ignores the s.coarsen() method...
	for (scalar s in listc)
	  s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
	for (vector v in listf)
	  foreach_dimension() {
	    v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
	    v.x[1,0] = (fine(v.x,2,0) + fine(v.x,2,1))/2.;
	  }
      }
      boundary_iterate (restriction, list, l);
    }
  }
  free (listc);
  free (listf);
}

void wavelet (scalar s, scalar w)
{
  restriction ({s});
  foreach_fine_to_coarse() {
    double sc[4];
    int c = 0;
    foreach_child()
      sc[c++] = s[];
    s.prolongation (point, s);
    c = 0;
    foreach_child() {
      /* difference between fine value and its prolongation */
      w[] = sc[c] - s[];
      s[] = sc[c++];
    }
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

  assert (fabs(4.*cm[] 
	       - fine(cm,0,0) - fine(cm,1,0) 
	       - fine(cm,0,1) - fine(cm,1,1)) < 1e-10);

  /* for each child */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      /* linear interpolation from coarser level (conservative) */
      fine(s,k,l) = s[] + 
	((2*k-1)*g.x*(fine(cm,!k,0) + fine(cm,!k,1)) +
	 (2*l-1)*g.y*(fine(cm,0,!l) + fine(cm,1,!l)))/(8.*cm[]);
}

void refine_reset (Point point, scalar v)
{
  /* foreach_child() */
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      fine(v,k,l) = 0.;
}

void refine_injection (Point point, scalar v)
{
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      fine(v,k,l) = v[];
}

void coarsen_face (Point point, scalar s)
{
  vector v = s.v;
  foreach_dimension() {
    v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
    v.x[1,0] = (fine(v.x,2,0) + fine(v.x,2,1))/2.;
  }
}

void no_coarsen (Point point, scalar s) {}

vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  v.x.coarsen = coarsen_face;
  v.y.coarsen = no_coarsen;
  return v;
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
    for (int k = -2; k <= 3; k++)
      for (int l = -2; l <= 3; l++) {
	for (scalar v in all) {
	  fprintf (fp, "%g %g ", 
		   xf + k*Delta/2. + v.d.x*Delta/4., 
		   yf + l*Delta/2. + v.d.y*Delta/4.);
	  if (allocated_child(k,l))
	    fprintf (fp, "%g ", fine(v,k,l));
	  else
	    fputs ("n/a ", fp);
	}
	fputc ('\n', fp);
      }
    fclose (fp);
    fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug                = multigrid_debug;
  init_face_vector     = multigrid_init_face_vector;
}
