#include "events.h"

void (* debug)    (Point);

@undef VARIABLES
@define VARIABLES							\
  double delta = L0*DELTA; /* cell size */				\
  /* cell/face center coordinates */					\
  double x  = (ig/2. + I + 0.5)*delta + X0;				\
  double y  = (jg/2. + J + 0.5)*delta + Y0;				\
  /* we need this to avoid compiler warnings */	                        \
  NOT_UNUSED(delta); NOT_UNUSED(x); NOT_UNUSED(y);			\
  /* and this when catching FPEs */					\
  _CATCH

#include "fpe.h"

@ifndef is_face_x
@ define is_face_x() true
@ define is_face_y() true
@endif

@ifndef foreach_boundary_ghost
@ define foreach_boundary_ghost(dir)					\
  foreach_boundary(dir,false) {						\
    point.i += ig; point.j += jg;					\
    ig = -ig; jg = -jg;							\
    POINT_VARIABLES;
@ define end_foreach_boundary_ghost()					\
    ig = -ig; jg = -jg;							\
  } end_foreach_boundary()
@endif

@define end_foreach_face()

#define boundary_flux(...)
#define output_stencil(v,fp) _output_stencil(point,v,#v,fp)
void _output_stencil (Point point, scalar s, const char * name, FILE * fp)
{
  int width = 25, len = (width - strlen(name))/2. - 1;
  for (int i = 0; i < len; i++) fputc ('-', fp);
  fputc (' ', fp); fputs (name, fp); fputc (' ', fp);
  for (int i = 0; i < len; i++) fputc ('-', fp);
  fputc ('\n', fp);
  for (int j = GHOSTS; j >= -GHOSTS; j--) {
    if (J + j >= - GHOSTS && J + j < _n + GHOSTS) {
      for (int i = - GHOSTS; i <= GHOSTS; i++)
	if (I + i >= - GHOSTS && I + i < _n + GHOSTS) {
	  fprintf (fp, "%5.g", s[i,j]);
	  if ((I + i < 0 || I + i >= _n) &&
	      (J + j < 0 || J + j >= _n))
	    fputs (":C ", fp);
	  else if (I + i < 0 || I + i >= _n ||
		   J + j < 0 || J + j >= _n)
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

void (* boundary) (scalar *);

void cartesian_boundary (scalar * list)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary (b, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = s.boundary[b] (point, s);
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
    delta /= 2.;
    fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	     x - delta, y - delta,
	     x - delta, y + delta,
	     x + delta, y + delta,
	     x + delta, y - delta,
	     x - delta, y - delta);
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
      fprintf (fp, "%g %g", x + k*delta, y + l*delta);
      for (scalar v in all)
	fprintf (fp, " %g", v[k,l]);
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
  init_scalar = cartesian_init_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  boundary    = cartesian_boundary;
  debug       = cartesian_debug;
}
