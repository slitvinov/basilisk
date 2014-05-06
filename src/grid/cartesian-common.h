#include "events.h"

void (* debug)    (Point);

@define _val_constant(a,k,l) ((const double) _constant[a -_NVARMAX])

@undef VARIABLES
@def VARIABLES
  double Delta = L0*DELTA; /* cell size */
  /* cell/face center coordinates */
  double x  = (ig/2. + I + 0.5)*Delta + X0;
  double y  = (jg/2. + J + 0.5)*Delta + Y0;
  /* we need this to avoid compiler warnings */
  NOT_UNUSED(Delta); NOT_UNUSED(x); NOT_UNUSED(y);
  /* and this when catching FPEs */
  _CATCH;
@

#include "fpe.h"

@ifndef is_face_x
@ define is_face_x() true
@ define is_face_y() true
@endif

@ifndef foreach_boundary_ghost
@ def foreach_boundary_ghost(dir)
  foreach_boundary(dir,false) {
    point.i += ig; point.j += jg;
    ig = -ig; jg = -jg;
    POINT_VARIABLES;
@
@ def end_foreach_boundary_ghost()
    ig = -ig; jg = -jg;
  } end_foreach_boundary()
@
@endif

@ifndef foreach_boundary_face_ghost
@ def foreach_boundary_face_ghost(dir)
  foreach_boundary_face(dir) {
    point.i += ig; point.j += jg;
    ig = -ig; jg = -jg;
    POINT_VARIABLES;
@
@ def end_foreach_boundary_face_ghost()
    ig = -ig; jg = -jg;
  } end_foreach_boundary_face()
@
@endif

@define end_foreach_face()
@define end_foreach_vertex()

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
  _method = realloc (_method, nvar*sizeof (Methods));
  _method[nvar - 1] = (Methods){{NULL, NULL, NULL}};
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

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  for (scalar s in l) {
    scalar c = new scalar;
    _method[c] = _method[s];
    c.name = NULL;
    list = list_append (list, c);
  }
  return list;
}

void delete (scalar * list)
{
  if (all == NULL) // everything has already been freed
    return;

  if (list == all) {
    for (scalar f in all) {
      free (f.name); f.name = NULL;
    }
    all[0] = -1;
    return;
  }

  trash (list);
  for (scalar f in list) {
    free (f.name); f.name = NULL;
    scalar * s = all;
    for (; *s >= 0 && *s != f; s++);
    if (*s == f)
      for (; *s >= 0; s++)
	s[0] = s[1];
  }
}

// Cartesian methods

void (* boundary_level)   (scalar *, int l);
void (* boundary_normal)  (vector *);
void (* boundary_tangent) (vector *);

void boundary_face (vector * list)
{
  foreach_dimension() {
    foreach_boundary (right,false)
      for (vector u in list)
	u.x[ghost] = u.x.boundary[right] (point, u.x);
    foreach_boundary (left,false)
      for (vector u in list)
	u.x[] = u.x.boundary[left] (point, u.x);
  }
  boundary_normal (list);
  boundary_tangent (list);
}

void boundary (scalar * list)
{
  scalar * listc = NULL;
  vector * lists = NULL;
  for (scalar s in list) {
    if (s.face)
      lists = vectors_add (lists, s.v);
    else
      listc = list_add (listc, s);
  }
  if (listc) {
    boundary_level (listc, -1);
    free (listc);
  }
  if (lists) {
    boundary_face (lists);
    free (lists);
  }
}

void cartesian_boundary_level (scalar * list, int l)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary (b, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = s.boundary[b] (point, s);
}

void cartesian_boundary_normal (vector * list)
{
  // nothing to do
}

void cartesian_boundary_tangent (vector * list)
{
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);
  foreach_dimension() {
    for (int d = top; d <= bottom; d++)
      foreach_boundary_face (d)
	for (vector v in listv)
	  v.x[ghost] = v.x.boundary[d] (point, v.x);
  }
  free (listv);
}

static double symmetry (Point point, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, scalar s)
{
  return -s[];
}

static double noflux (Point point, scalar s)
{
  return 0.;
}

scalar cartesian_init_scalar (scalar s, const char * name)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    s.boundary[b] = s.boundary_homogeneous[b] = symmetry;
  s.prolongation = NULL; 
  s.refine = s.coarsen = NULL;
  s.gradient = NULL;
  if (name)
    s.name = strdup (name);
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
    v.x.boundary[right] = v.x.boundary[left] = noflux;  
    v.x.d.x = -1;
    v.x.face = true;
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
  for (int k = -1; k <= 1; k++)
    for (int l = -1; l <= 1; l++) {
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
  init_scalar       = cartesian_init_scalar;
  init_vector       = cartesian_init_vector;
  init_tensor       = cartesian_init_tensor;
  boundary_level    = cartesian_boundary_level;
  boundary_normal   = cartesian_boundary_normal;
  boundary_tangent  = cartesian_boundary_tangent;
  debug             = cartesian_debug;
  init_face_vector  = cartesian_init_face_vector;
}
