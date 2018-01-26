#include "events.h"

void (* debug)    (Point);

@define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

@undef VARIABLES
@def VARIABLES
  double Delta = L0*_DELTA; /* cell size */
  foreach_dimension()
    double Delta_x = Delta; /* cell size (with mapping) */
  /* cell/face center coordinates */
  double x = (ig/2. + _I + 0.5)*Delta + X0; NOT_UNUSED(x);
#if dimension > 1
  double y = (jg/2. + _J + 0.5)*Delta + Y0;
#else
  double y = 0.;
#endif
 NOT_UNUSED(y);
#if dimension > 2
  double z = (kg/2. + _K + 0.5)*Delta + Z0;
#else
  double z = 0.;
#endif
  NOT_UNUSED(z);
  /* we need this to avoid compiler warnings */
  NOT_UNUSED(Delta);
  foreach_dimension()
    NOT_UNUSED(Delta_x);
  /* and this when catching FPEs */
  _CATCH;
@

#include "fpe.h"

@define end_foreach_face()

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  scalar s;
  for (s.i = 0; s.i < nvar; s.i++)
    if (!list_lookup (all, s)) { // found a previously freed slot
      init_scalar (s, name);
      trash (((scalar []){s, {-1}}));
      all = list_append (all, s);
      return s;
    }
  
  // need to allocate a new slot
  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  qrealloc (_attribute, nvar, _Attributes);
  memset (&_attribute[nvar-1], 0, sizeof (_Attributes));
  s = (scalar){nvar - 1};
  realloc_scalar(); // allocate extra space on the grid
  init_scalar (s, name);
  trash (((scalar []){s, {-1}}));
  all = list_append (all, s);
  return s;
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_scalar (name);
  foreach_dimension()
    s.d.x = -1;
  return s;
}

static vector alloc_vector (const char * name)
{
  vector v;
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  foreach_dimension() {
    sprintf (cname, ext.x, name);
    v.x = new_scalar (cname);
  }
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_face_vector (v, NULL);
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  foreach_dimension() {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  foreach_dimension() {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
  #if dimension > 1
    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
  #endif
  #if dimension > 2
    sprintf (cname, "%s.x.z", name);
    t.x.z = new_scalar(cname);
    t.z.x = t.x.z;
    sprintf (cname, "%s.y.z", name);
    t.y.z = new_scalar(cname);
    t.z.y = t.y.z;
  #endif
  /* fixme: boundary conditions don't work!  This is because boundary
     attributes depend on the index and should (but cannot) be
     different for t.x.y and t.y.x */
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    qrealloc (_constant, nconst, double);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  foreach_dimension()
    init_const_scalar (v.x, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  foreach_dimension()
    v.x.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = a.name;
  double (** boundary) (Point, Point, scalar) = a.boundary;
  double (** boundary_homogeneous) (Point, Point, scalar) =
    a.boundary_homogeneous;
  _attribute[a.i] = _attribute[b.i];
  a.name = name;
  a.boundary = boundary;
  a.boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    a.boundary[i] = b.boundary[i];
    a.boundary_homogeneous[i] = b.boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  for (scalar s in l) {
    scalar c = new scalar;
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  for (scalar s in list)
    foreach_dimension()
      if (s.v.x.i >= 0 && map[s.v.x.i] >= 0)
	s.v.x.i = map[s.v.x.i];
  return list;
}

void delete (scalar * list)
{
  if (all == NULL) // everything has already been freed
    return;

  for (scalar f in list) {
    if (f.delete)
      f.delete (f);
    free (f.name); f.name = NULL;
    free (f.boundary); f.boundary = NULL;
    free (f.boundary_homogeneous); f.boundary_homogeneous = NULL;
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  for (scalar f in list) {
    scalar * s = all;
    for (; s->i >= 0 && s->i != f.i; s++);
    if (s->i == f.i)
      for (; s->i >= 0; s++)
	s[0] = s[1];
  }
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

void free_solver()
{
  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }
  
  delete (all);
  free (all); all = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      free (e);
      e = next;
    }
  }

  free (Events); Events = NULL;
  free (_attribute); _attribute = NULL;
  free (_constant); _constant = NULL;
  free_grid();
  qpclose_all();
@if TRACE
  trace_off();
@endif
@if MTRACE
  pmuntrace();
@endif
@if _CADNA
  cadna_end();
@endif
}

// Cartesian methods

void (* boundary_level) (scalar *, int l);
void (* boundary_flux)  (vector *);

trace
void boundary (scalar * list)
{
  if (list == NULL)
    return;
  vector * listf = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.face)
      listf = vectors_add (listf, s.v);
  if (listf) {
    boundary_flux (listf);
    free (listf);
  }
  boundary_level (list, -1);
}

void cartesian_boundary_level (scalar * list, int l)
{
  boundary_iterate (level, list, l);
}

void cartesian_boundary_flux (vector * list)
{
  // nothing to do
}

static double symmetry (Point point, Point neighbor, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, Point neighbor, scalar s)
{
  return -s[];
}

double (* default_scalar_bc[]) (Point, Point, scalar) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{
  // keep name
  char * pname;
  if (name) {
    free (s.name);
    pname = strdup (name);
  }
  else
    pname = s.name;
  free (s.boundary);
  free (s.boundary_homogeneous);
  // reset all attributes
  _attribute[s.i] = (const _Attributes){0};
  s.name = pname;
  /* set default boundary conditions */
  s.boundary = (double (**)(Point, Point, scalar))
    malloc (nboundary*sizeof (void (*)()));
  s.boundary_homogeneous = (double (**)(Point, Point, scalar))
    malloc (nboundary*sizeof (void (*)()));
  for (int b = 0; b < nboundary; b++)
    s.boundary[b] = s.boundary_homogeneous[b] =
      b < 2*dimension ? default_scalar_bc[b] : symmetry;
  s.gradient = NULL;
  foreach_dimension() {
    s.d.x = 0;  // not face
    s.v.x.i = -1; // not a vector component
  }
  s.face = false;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  foreach_dimension() {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    v.x.v = v;
  }
  /* set default boundary conditions */
  for (int d = 0; d < nboundary; d++)
    v.x.boundary[d] = v.x.boundary_homogeneous[d] =
      d < 2*dimension ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  foreach_dimension() {
    v.x.d.x = -1;
    v.x.face = true;
  }
  for (int d = 0; d < nboundary; d++)
    v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  foreach_dimension() {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
  /* set default boundary conditions */
  #if dimension == 1
    for (int b = 0; b < nboundary; b++)
      t.x.x.boundary[b] = t.x.x.boundary_homogeneous[b] =
	b < 2*dimension ? default_scalar_bc[b] : symmetry;
  #elif dimension == 2
    for (int b = 0; b < nboundary; b++) {
      t.x.x.boundary[b] = t.y.x.boundary[b] = 
	t.x.x.boundary_homogeneous[b] = t.y.y.boundary_homogeneous[b] = 
	b < 2*dimension ? default_scalar_bc[b] : symmetry;
      t.x.y.boundary[b] = t.y.y.boundary[b] = 
	t.x.y.boundary_homogeneous[b] = t.y.x.boundary_homogeneous[b] = 
	b < 2*dimension ? default_vector_bc[b] : antisymmetry;
    }
  #else
    assert (false); // not implemented yet
  #endif
  return t;
}

void output_cells (FILE * fp)
{
  foreach() {
    Delta /= 2.;
    #if dimension == 1
      fprintf (fp, "%g 0\n%g 0\n\n", x - Delta, x + Delta);
    #elif dimension == 2
      fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	       x - Delta, y - Delta,
	       x - Delta, y + Delta,
	       x + Delta, y + Delta,
	       x + Delta, y - Delta,
	       x - Delta, y - Delta);
    #else // dimension == 3
      for (int i = -1; i <= 1; i += 2) {
	fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
		 x - Delta, y - Delta, z + i*Delta,
		 x - Delta, y + Delta, z + i*Delta,
		 x + Delta, y + Delta, z + i*Delta,
		 x + Delta, y - Delta, z + i*Delta,
		 x - Delta, y - Delta, z + i*Delta);
	for (int j = -1; j <= 1; j += 2)
	  fprintf (fp, "%g %g %g\n%g %g %g\n\n",
		   x + i*Delta, y + j*Delta, z - Delta,
		   x + i*Delta, y + j*Delta, z + Delta);
      }
    #endif
  }
  fflush (fp);
}

static char * replace_ (const char * vname)
{
  char * name = strdup (vname), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
			const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp, 
	   "  load 'debug.plot'\n"
	   "  v=%s\n"
#if dimension == 1   
	   "  plot '%s' w l lc 0, "
	   "'%s' u 1+2*v:(0):2+2*v w labels tc lt 1 title columnhead(2+2*v)",
#elif dimension == 2
	   "  plot '%s' w l lc 0, "
	   "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",
#elif dimension == 3
	   "  splot '%s' w l lc 0, "
	   "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",
#endif
	   vname, cells, stencil);
  free (vname);
}

void cartesian_debug (Point point)
{
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  for (scalar v in all)
#if dimension == 1
    fprintf (fp, "x %s ", v.name);
#elif dimension == 2
    fprintf (fp, "x y %s ", v.name);
#elif dimension == 3
    fprintf (fp, "x y z %s ", v.name);
#endif
  fputc ('\n', fp);
  #if dimension == 1
    for (int k = -2; k <= 2; k++) {
      for (scalar v in all) {
	fprintf (fp, "%g ", x + k*Delta + v.d.x*Delta/2.);
	if (allocated(k))
	  fprintf (fp, "%g ", v[k]);
	else
	  fputs ("n/a ", fp);
      }
      fputc ('\n', fp);
    }
  #elif dimension == 2
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
	for (scalar v in all) {
	  fprintf (fp, "%g %g ",
		   x + k*Delta + v.d.x*Delta/2., 
		   y + l*Delta + v.d.y*Delta/2.);
	  if (allocated(k,l))
	    fprintf (fp, "%g ", v[k,l]);
	  else
	    fputs ("n/a ", fp);
	}
	fputc ('\n', fp);
      }
  #elif dimension == 3
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
	for (int m = -2; m <= 2; m++) {
	  for (scalar v in all) {
	    fprintf (fp, "%g %g %g ",
		     x + k*Delta + v.d.x*Delta/2., 
		     y + l*Delta + v.d.y*Delta/2.,
		     z + m*Delta + v.d.z*Delta/2.);
	    if (allocated(k,l,m))
	      fprintf (fp, "%g ", v[k,l,m]);
	    else
	      fputs ("n/a ", fp);
	  }
	  fputc ('\n', fp);
	}
  #endif
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp, 
	   "set term x11\n"
	   "set size ratio -1\n"
	   "set key outside\n");
  for (scalar s in all) {
    char * name = replace_ (s.name);
    fprintf (fp, "%s = %d\n", name, s.i);
    free (name);
  }
  fclose (fp);

  fprintf (stderr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (stderr, _attribute[0].name, name, stencil);
  fflush (stderr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar      = cartesian_init_scalar;
  init_vector      = cartesian_init_vector;
  init_tensor      = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level   = cartesian_boundary_level;
  boundary_flux    = cartesian_boundary_flux;
  debug            = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{
  scalar v = p.v;
#if dimension == 1
  x = (p.x - x)/Delta - v.d.x/2.;
  int i = sign(x);
  x = fabs(x);
  /* linear interpolation */
  return v[]*(1. - x) + v[i]*x;
#elif dimension == 2
  x = (p.x - x)/Delta - v.d.x/2.;
  y = (p.y - y)/Delta - v.d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);
  /* bilinear interpolation */
  return ((v[]*(1. - x) + v[i]*x)*(1. - y) + 
	  (v[0,j]*(1. - x) + v[i,j]*x)*y);
#else // dimension == 3
  x = (p.x - x)/Delta - v.d.x/2.;
  y = (p.y - y)/Delta - v.d.y/2.;
  z = (p.z - z)/Delta - v.d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);
  /* trilinear interpolation */
  return (((v[]*(1. - x) + v[i]*x)*(1. - y) + 
	   (v[0,j]*(1. - x) + v[i,j]*x)*y)*(1. - z) +
	  ((v[0,0,k]*(1. - x) + v[i,0,k]*x)*(1. - y) + 
	   (v[0,j,k]*(1. - x) + v[i,j,k]*x)*y)*z);
#endif  
}

trace
double interpolate (struct _interpolate p)
{
  Point point = locate (p.x, p.y, p.z);
  if (point.level < 0)
    return nodata;
  return interpolate_linear (point, p);
}

trace
void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate (a[i].x, a[i].y, a[i].z);
    if (point.level >= 0) {
      for (scalar s in list)
	v[j++] = !linear ? s[] :
	  interpolate_linear (point,
			      (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      for (scalar s in list)
	v[j++] = nodata;
  }
@if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
@endif
}

// Boundaries

typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  for (scalar s in all) {
    s.boundary = (double (**)(Point, Point, scalar))
      realloc (s.boundary, nboundary*sizeof (void (*)()));
    s.boundary_homogeneous = (double (**)(Point, Point, scalar))
      realloc (s.boundary_homogeneous, nboundary*sizeof (void (*)()));
  }
  for (scalar s in all) {
    if (s.v.x.i < 0) // scalar
      s.boundary[b] = s.boundary_homogeneous[b] = symmetry;
    else if (s.v.x.i == s.i) { // vector
      vector v = s.v;
      foreach_dimension()
	v.y.boundary[b] = v.y.boundary_homogeneous[b] = symmetry;
      v.x.boundary[b] = v.x.boundary_homogeneous[b] =
	v.x.face ? NULL : antisymmetry;
    }
  }
  return b;
}

// Periodic boundary conditions

static double periodic_bc (Point point, Point neighbor, scalar s)
{
  return nodata; // should not be used
}

static void periodic_boundary (int d)
{
  /* We change the conditions for existing scalars. */
  for (scalar s in all)
    s.boundary[d] = s.boundary_homogeneous[d] = periodic_bc;
  /* Normal components of face vector fields should remain NULL. */
  for (scalar s in all)
    if (s.face) {
      vector v = s.v;
      v.x.boundary[d] = v.x.boundary_homogeneous[d] = NULL;
    }
  /* We also change the default boundary conditions (for new fields). */
  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

int Periods[dimension];

void periodic (int dir)
{
  #if dimension < 2
    assert (dir <= left);
  #elif dimension < 3
    assert (dir <= bottom);
  #else
    assert (dir <= back);
  #endif
  // This is the component in the given direction i.e. 0 for x and 1 for y
  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  Periods[c] = true;
}
