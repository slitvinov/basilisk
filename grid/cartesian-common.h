#include "events.h"

void (* debug)    (Point);

// coordinates of the center of the box
double X0 = -0.5, Y0 = -0.5;
// size of the box
double L0 = 1.;
#define DX    (L0*delta)
#define XC(i) ((i + 0.5)*DX + X0*L0)
#define YC(j) ((j + 0.5)*DX + X0*L0)

#undef VARIABLES
#define VARIABLES							\
  double delta = DELTA;          /* cell size */			\
  double x  = XC(ig/2. + I), y  = YC(jg/2. + J); /* cell/face center */	\
  /* we need this to avoid compiler warnings */	                        \
  NOT_UNUSED(delta); NOT_UNUSED(x); NOT_UNUSED(y);			\
  /* and this when catching FPEs */					\
  _CATCH

#include "fpe.h"

#ifndef is_face_x
# define is_face_x() true
# define is_face_y() true
#endif

#ifndef foreach_boundary_ghost
# define foreach_boundary_ghost(dir)					\
  foreach_boundary(dir,false) {						\
    point.i += ig; point.j += jg;					\
    ig = -ig; jg = -jg;							\
    POINT_VARIABLES;
# define end_foreach_boundary_ghost()					\
    ig = -ig; jg = -jg;							\
  } end_foreach_boundary()
#endif

#define boundary_ghost(d, x) {						\
    foreach_boundary_ghost (d) { x; } end_foreach_boundary_ghost();	\
  }

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

// Cartesian methods

void (* boundary) (scalar *);

void cartesian_boundary (scalar * list)
{
  for (int b = 0; b < nboundary; b++)
    foreach_boundary (b, true) // also traverse corners
      for (scalar s in list)
	s[ghost] = _boundary[b][s] (point, s);
}

static double symmetry (Point point, scalar s)
{
  return s[];
}

static double antisymmetry (Point point, scalar s)
{
  return -s[];
}

scalar cartesian_new_scalar (scalar s)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++)
    _boundary[b][s] = symmetry;
  return s;
}

vector cartesian_new_vector (vector v)
{
  /* set default boundary conditions (symmetry) */
  _boundary[top][v.x] = _boundary[bottom][v.x] = symmetry;
  _boundary[right][v.y] = _boundary[left][v.y] = symmetry;
  _boundary[right][v.x] = _boundary[left][v.x] = antisymmetry;
  _boundary[top][v.y] = _boundary[bottom][v.y] = antisymmetry;
  return v;
}

tensor cartesian_new_tensor (tensor t)
{
  /* set default boundary conditions (symmetry) */
  for (int b = 0; b < nboundary; b++) {
    _boundary[b][t.x.x] = _boundary[b][t.y.y] = symmetry;
    _boundary[b][t.x.y] = _boundary[b][t.y.x] = antisymmetry;
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
      for (scalar v = 0; v < nvar; v++)
	fprintf (fp, " %g", val(v,k,l));
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
  new_scalar = cartesian_new_scalar;
  new_vector = cartesian_new_vector;
  new_tensor = cartesian_new_tensor;
  boundary   = cartesian_boundary;
  debug = cartesian_debug;
}
