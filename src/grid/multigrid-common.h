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

static inline void coarsen_average (Point point, scalar s)
{
  double sum = 0.;
  foreach_child()
    sum += s[];
  s[] = sum/(1 << dimension);
}

static inline void coarsen_volume_average (Point point, scalar s)
{
  double sum = 0.;
  foreach_child()
    sum += cm[]*s[];
  s[] = sum/(1 << dimension)/cm[];
}

static inline void face_average (Point point, vector v)
{
  foreach_dimension() {
    #if dimension == 1
      v.x[] = fine(v.x,0);
      v.x[1] = fine(v.x,2);
    #elif dimension == 2
      v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
      v.x[1] = (fine(v.x,2,0) + fine(v.x,2,1))/2.;
    #else // dimension == 3
      v.x[] = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
	       fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      v.x[1] = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
		fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;
    #endif
  }
}
  
static inline void coarsen_face (Point point, scalar s)
{
  face_average (point, s.v);
}

static inline void no_coarsen (Point point, scalar s) {}

static inline void no_data (Point point, scalar s) {
  foreach_child()
    s[] = nodata;
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
	  coarsen_average (point, s);
	for (vector v in listf)
	  face_average (point, v);
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
    double sc[1 << dimension];
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

static inline double bilinear (Point point, scalar s)
{
  #if dimension == 1
    return (3.*coarse(s,0) + coarse(s,child.x))/4.;
  #elif dimension == 2
    return (9.*coarse(s,0) + 
	    3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
	    coarse(s,child.x,child.y))/16.;
  #else // dimension == 3
    return (27.*coarse(s,0) + 
	    9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		coarse(s,0,0,child.z)) + 
	    3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		coarse(s,0,child.y,child.z)) + 
	    coarse(s,child.x,child.y,child.z))/64.;
  #endif
}

static inline void refine_bilinear (Point point, scalar s)
{
  foreach_child()
    s[] = bilinear (point, s);
}

static inline double biquadratic (Point point, scalar s)
{
#if dimension == 1
  assert (false);
  return 0.;
#elif dimension == 2
  return (900.*coarse(s,0,0) + 25.*coarse(s,child.x,child.y) +
	  150.*(coarse(s,child.x,0) + coarse(s,0,child.y)) +
	  9.*coarse(s,-child.x,-child.y) -
	  90.*(coarse(s,-child.x, 0) + coarse(s,0,-child.y)) -
	  15.*(coarse(s,child.x,-child.y) + coarse(s,-child.x,child.y)))/1024.;
#else // dimension == 3
  assert (false);
  return 0.;
#endif
}

static inline void refine_biquadratic (Point point, scalar s)
{
  foreach_child()
    s[] = biquadratic (point, s);
}

static inline void refine_linear (Point point, scalar s)
{
  coord g;
  if (s.gradient)
    foreach_dimension()
      g.x = s.gradient (s[-1], s[], s[1]);
  else
    foreach_dimension()
      g.x = (s[1] - s[-1])/2.;

  double sc = s[], cmc = 4.*cm[], sum = cm[]*(1 << dimension);
  foreach_child() {
    s[] = sc;
    foreach_dimension()
      s[] += child.x*g.x*cm[-child.x]/cmc;
    sum -= cm[];
  }
  assert (fabs(sum) < 1e-10);
}

static inline void refine_reset (Point point, scalar v)
{
  foreach_child()
    v[] = 0.;
}

static inline void refine_injection (Point point, scalar v)
{
  double val = v[];
  foreach_child()
    v[] = val;
}

vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  foreach_dimension()
    v.y.coarsen = no_coarsen;
  v.x.coarsen = coarsen_face;
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
    #if dimension == 1
      double xc = x - child.x*Delta/2.;
      for (int k = 0; k <= 1; k++)
	for (scalar v in all)
	  fprintf (fp, "%g %g ", 
		   xc + k*child.x*Delta*2. + v.d.x*Delta, 
		   coarse(v,k*child.x));
      fputc ('\n', fp);
      fprintf (stderr, ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 3 t ''", name);
    #elif dimension == 2
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
      fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
    #elif dimension == 3
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
	for (int l = 0; l <= 1; l++)
	  for (int m = 0; m <= 1; m++) {
	    for (scalar v in all)
	      fprintf (fp, "%g %g %g %g ", 
		       xc + k*child.x*Delta*2. + v.d.x*Delta, 
		       yc + l*child.y*Delta*2. + v.d.y*Delta,
		       zc + m*child.z*Delta*2. + v.d.z*Delta,
		       coarse(v,k*child.x,l*child.y,m*child.z));
	    fputc ('\n', fp);
	  }
      fprintf (stderr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
	       name);
    #endif
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
    #if dimension == 1
      double xf = x - Delta/4.;
      for (int k = -2; k <= 3; k++)
	for (scalar v in all) {
	  fprintf (fp, "%g ", xf + k*Delta/2. + v.d.x*Delta/4.);
	  if (allocated_child(k))
	    fprintf (fp, "%g ", fine(v,k));
	  else
	    fputs ("n/a ", fp);
	}
      fputc ('\n', fp);
      fprintf (stderr, ", '%s' u 1+2*v:(0):2+2*v w labels tc lt 2 t ''", name);
    #elif dimension == 2
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
      fprintf (stderr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
    #elif dimension == 3
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
	for (int l = -2; l <= 3; l++)
	  for (int m = -2; m <= 3; m++) {
	    for (scalar v in all) {
	      fprintf (fp, "%g %g %g ", 
		       xf + k*Delta/2. + v.d.x*Delta/4., 
		       yf + l*Delta/2. + v.d.y*Delta/4.,
		       zf + m*Delta/2. + v.d.z*Delta/4.);
	      if (allocated_child(k,l,m))
		fprintf (fp, "%g ", fine(v,k,l,m));
	      else
		fputs ("n/a ", fp);
	    }
	    fputc ('\n', fp);
	  }
      fprintf (stderr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
	       name);
    #endif
    fclose (fp);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug                = multigrid_debug;
  init_face_vector     = multigrid_init_face_vector;
}
