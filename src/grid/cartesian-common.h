#include "events.h"

void (* debug)    (Point);

@define _val_constant(a,k,l) ((const double) _constant[a -_NVARMAX])

@undef VARIABLES
@def VARIABLES
  double Delta = L0*DELTA; /* cell size */
  double Delta_x = Delta, Delta_y = Delta; /* cell size (with mapping) */
  /* cell/face center coordinates */
  double x = (ig/2. + I + 0.5)*Delta + X0;
  double y = (jg/2. + J + 0.5)*Delta + Y0;
  /* we need this to avoid compiler warnings */
  NOT_UNUSED(Delta); NOT_UNUSED(Delta_x); NOT_UNUSED(Delta_y);
  NOT_UNUSED(x); NOT_UNUSED(y);
  /* and this when catching FPEs */
  _CATCH;
@

#include "fpe.h"

@define end_foreach_face()

#define output_stencil(v,fp) _output_stencil(point,v,#v,fp)
void _output_stencil (Point point, scalar s, const char * name, FILE * fp)
{
  int width = 25, len = (width - strlen(name))/2. - 1;
  for (int i = 0; i < len; i++) fputc ('-', fp);
  fputc (' ', fp); fputs (name, fp); fputc (' ', fp);
  for (int i = 0; i < len; i++) fputc ('-', fp);
  fputc ('\n', fp);
  for (int j = GHOSTS; j >= -GHOSTS; j--) {
    if (J + j >= - GHOSTS && J + j < NN + GHOSTS) {
      for (int i = - GHOSTS; i <= GHOSTS; i++)
	if (I + i >= - GHOSTS && I + i < NN + GHOSTS) {
	  fprintf (fp, "%5.g", s[i,j]);
	  if ((I + i < 0 || I + i >= NN) &&
	      (J + j < 0 || J + j >= NN))
	    fputs (":C ", fp);
	  else if (I + i < 0 || I + i >= NN ||
		   J + j < 0 || J + j >= NN)
	    fputs (":B ", fp);
	  else
	    fputs ("   ", fp);
	}
	else
	  fputs ("   ?    ", fp);
    }
    else
      fputs ("???????????????????????", fp);
    fputc ('\n', fp);
  }
}

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  for (int i = 0; i < nvar; i++)
    if (!list_lookup (all, i)) { // found a previously freed slot
      all = list_append (all, i);
      init_scalar (i, name);
      trash (((scalar []){i, -1}));
      return i;
    }
  
  // need to allocate a new slot
  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  _attribute = realloc (_attribute, nvar*sizeof (_Attributes));
  all = realloc (all, sizeof (scalar)*(nvar + 1));
  all[nvar - 1] = nvar - 1;
  all[nvar] = -1;
  realloc_scalar(); // allocate extra space on the grid
  init_scalar (nvar - 1, name);
  trash (((scalar []){nvar - 1, -1}));
  return nvar - 1;
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
  sprintf (cname, "%s.x", name);
  v.x = new_scalar (cname);
  sprintf (cname, "%s.y", name);
  v.y = new_scalar (cname);
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
  sprintf (cname, "%s.x", name);
  vector tx = new_vector (cname);
  sprintf (cname, "%s.y", name);
  vector ty = new_vector (cname);
  tensor t;
  t.x = tx; t.y = ty;
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  tensor t;
  sprintf (cname, "%s.x.x", name);
  t.x.x = new_scalar(cname);
  sprintf (cname, "%s.x.y", name);
  t.x.y = new_scalar(cname);
  sprintf (cname, "%s.y.y", name);
  t.y.y = new_scalar(cname);
  t.y.x = t.x.y;
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s - _NVARMAX >= nconst) {
    nconst = s - _NVARMAX + 1;
    _constant = realloc (_constant, nconst*sizeof (double));
  }
  _constant[s - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = i + _NVARMAX;
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  init_const_scalar (v.x, name, val[0]);
  init_const_scalar (v.y, name, val[1]);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v = {i + _NVARMAX, i + 1 + _NVARMAX};
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  free (a.name);
  _attribute[a] = _attribute[b];
  a.name = NULL;
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  for (scalar s in l) {
    scalar c = new scalar;
    scalar_clone (c, s);
    list = list_append (list, c);
  }
  return list;
}

void delete (scalar * list)
{
  if (all == NULL) // everything has already been freed
    return;

  for (scalar f in list) {
    free (f.name); f.name = NULL;
  }

  if (list == all) {
    all[0] = -1;
    return;
  }

  trash (list);
  for (scalar f in list) {
    scalar * s = all;
    for (; *s >= 0 && *s != f; s++);
    if (*s == f)
      for (; *s >= 0; s++)
	s[0] = s[1];
  }
}

void free_solver()
{
  delete (all);
  free (all); all = NULL;
  free (Events); Events = NULL;
  free (_attribute); _attribute = NULL;
  free_grid();
@if TRACE
  trace_off();
@endif
}

// Cartesian methods

void (* boundary_level) (scalar *, int l);
void (* boundary_flux)  (vector *);

void boundary (scalar * list)
{
  vector * listf = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.face)
      listf = vectors_add (listf, s.v);
  if (listf)
    boundary_flux (listf);
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

static double symmetry (Point point, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, scalar s)
{
  return -s[];
}

scalar cartesian_init_scalar (scalar s, const char * name)
{
  // keep name
  char * pname = name ? strdup (name) : s.name;
  // reset all attributes
  _attribute[s] = (const _Attributes){{0}};
  s.name = pname;
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    s.boundary[b] = s.boundary_homogeneous[b] = symmetry;
  s.gradient = NULL;
  foreach_dimension() {
    s.d.x = 0;  // not face
    s.v.x = -1; // not a vector component
  }
  s.face = false;
  return s;
}

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
    /* set default boundary conditions (symmetry) */
    v.x.boundary[right] = v.x.boundary[left] = 
      v.x.boundary_homogeneous[right] = v.x.boundary_homogeneous[left] = 
      antisymmetry;
    v.x.v = v;
  }
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  foreach_dimension() {
    v.x.boundary[right] = v.x.boundary[left] = NULL;
    v.x.d.x = -1;
    v.x.face = v.x.normal = true;
  }
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
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    t.x.x.boundary[b] = t.y.y.boundary[b] = 
      t.x.x.boundary_homogeneous[b] = t.y.y.boundary_homogeneous[b] = 
      symmetry;
    t.x.y.boundary[b] = t.y.x.boundary[b] = 
      t.x.y.boundary_homogeneous[b] = t.y.x.boundary_homogeneous[b] = 
      antisymmetry;
  }
  return t;
}

void output_cells (FILE * fp)
{
  foreach() {
    Delta /= 2.;
    fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	     x - Delta, y - Delta,
	     x - Delta, y + Delta,
	     x + Delta, y + Delta,
	     x + Delta, y - Delta,
	     x - Delta, y - Delta);
  }
  fflush (fp);
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
    fprintf (fp, "x y %s ", v.name);
  fputc ('\n', fp);
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
  fclose (fp);

  fprintf (stderr, 
	   "Last point stencils can be displayed using (in gnuplot)\n"
	   "  set size ratio -1\n"
	   "  set key outside\n"
	   "  v=0\n"
	   "  plot '%s' w l lc 0, "
	   "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",
	   name, stencil);
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

Point locate (double xp, double yp);

double interpolate (scalar v, double xp, double yp)
{
  Point point = locate (xp, yp);
  if (point.level < 0)
    return nodata;
  x = (xp - x)/Delta - v.d.x/2.;
  y = (yp - y)/Delta - v.d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);
  /* bilinear interpolation */
  return (v[]*(1. - x)*(1. - y) + 
	  v[i,0]*x*(1. - y) + 
	  v[0,j]*(1. - x)*y + 
	  v[i,j]*x*y);
}
