#include "events.h"

void (* debug)    (Point);

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
      return i;
    }
  
  // need to allocate a new slot
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

vector new_vector (const char * name)
{
  char cname[strlen(name) + 3];
  sprintf (cname, "%s.x", name);
  scalar vx = new_scalar (cname);
  sprintf (cname, "%s.y", name);
  scalar vy = new_scalar (cname);
  vector v;
  v.x = vx; v.y = vy;
  init_vector (v, name);
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
  init_tensor (t, name);
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
  init_tensor (t, name);
  return t;
}

scalar * clone (scalar * l)
{
  scalar * list = NULL;
  for (scalar s in l) {
    scalar c = new scalar;
    _method[c] = _method[s];
    list = list_append (list, c);
  }
  return list;
}

void delete (scalar * list)
{
  if (all == NULL) // everything has already been freed
    return;
  trash (list);
  for (scalar f in list) {
    scalar * s = all;
    for (; *s >= 0 && *s != f; s++);
    if (*s == f)
      for (; *s >= 0; s++)
	s[0] = s[1];
  }
}

// Cartesian methods

void (* boundary)         (scalar *);
void (* boundary_normal)  (vector *);
void (* boundary_tangent) (vector *);

void cartesian_boundary (scalar * list)
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
  // fixme: hardcoded for 2D
  for (int d = top; d <= bottom; d++)
    foreach_boundary_face (d)
      for (vector v in list)
	v.x[ghost] = v.x.boundary[d] (point, v.x);
  for (int d = right; d <= left; d++)
    foreach_boundary_face (d)
      for (vector v in list)
	v.y[ghost] = v.y.boundary[d] (point, v.y);
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
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    s.boundary[b] = symmetry;
  return s;
}

vector cartesian_init_vector (vector v, const char * name)
{
  foreach_dimension()
    init_scalar (v.x, name);
  /* set default boundary conditions (symmetry) */
  v.x.boundary[right] = v.x.boundary[left] = antisymmetry;
  v.y.boundary[top] = v.y.boundary[bottom] = antisymmetry;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  foreach_dimension()
    init_vector (t.x, name);
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    t.x.x.boundary[b] = t.y.y.boundary[b] = symmetry;
    t.x.y.boundary[b] = t.y.x.boundary[b] = antisymmetry;
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
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);

  fp = fopen ("stencil", "w");
  for (int k = -1; k <= 1; k++)
    for (int l = -1; l <= 1; l++) {
      fprintf (fp, "%g %g", x + k*Delta, y + l*Delta);
      if (allocated (k,l)) {
	for (scalar v in all)
	  fprintf (fp, " %g", v[k,l]);
      }
      else
	for (scalar v in all)
	  fprintf (fp, " n/a");
      fputc ('\n', fp);
    }
  fclose (fp);

  fputs ("Last point stencils can be displayed using e.g.\n"
	 "gnuplot> v=0\n"
	 "gnuplot> plot 'cells' w l lc 0, "
	 "'stencil' u 1:2:3+v w labels tc lt 1",
	 stderr);
}

void cartesian_methods()
{
  init_scalar      = cartesian_init_scalar;
  init_vector      = cartesian_init_vector;
  init_tensor      = cartesian_init_tensor;
  boundary         = cartesian_boundary;
  boundary_normal  = cartesian_boundary_normal;
  boundary_tangent = cartesian_boundary_tangent;
  debug            = cartesian_debug;
}

void boundary_mac (vector * list)
{
  foreach_boundary (right,false)
    for (vector u in list)
      u.x[ghost] = u.x.boundary[right] (point, u.x);
  foreach_boundary (top,false)
    for (vector u in list)
      u.y[ghost] = u.y.boundary[top] (point, u.y);
  foreach_boundary (left,false)
    for (vector u in list)
      u.x[] = u.x.boundary[left] (point, u.x);
  foreach_boundary (bottom,false)
    for (vector u in list)
      u.y[] = u.y.boundary[bottom] (point, u.y);
  boundary_normal (list);
  boundary_tangent (list);
}
