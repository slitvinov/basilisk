/**
# Automatic stencils and boundary conditions

Basilisk automatically computes, at runtime, the access pattern
(i.e. "stencils") of (basic) foreach loops (foreach(), foreach_face(),
foreach_vertex()).

This is done in practice by `qcc` which automatically adds, before
each foreach loop, a modified version of the loop body applied to a
single grid point of a special stencil grid.

The resulting access pattern is stored in the `read` and `write`
arrays associated with each field.

The `dirty` attribute is used to store the status of boundary
conditions for each field. */

attribute {
  double * write;
  int * read;
  int dirty; // // boundary conditions status:
  // 0: all conditions applied
  // 1: nothing applied
  // 2: boundary_face applied
}

/**
By default the stencil size is $5^\text{dimension}$.  */

#define _STENCIL 5
#if dimension == 1
# define _STENCIL_SIZE _STENCIL
#elif dimension == 2
# define _STENCIL_SIZE (_STENCIL*_STENCIL)
#elif dimension == 3
# define _STENCIL_SIZE (_STENCIL*_STENCIL*_STENCIL)
#endif

/**
This data structure contains information about the foreach loop being
processed. */

typedef struct {
  const char * fname; // name of the source file
  int line;           // line number in the source
  int first;          // is this the first time the loop is called?
  const char * each;  // the name of the foreach loop (foreach, foreach_face...)
  int face;           // the face component(s) being traversed
  bool vertex;        // is this a vertex traversal?
} ForeachData;

/**
## Stencil access detection

Given a (1D, 2D, 3D) index `i` this function returns the corresponding
index in the (1D) attribute arrays `read` and `write`, or -1 in case
of stencil overflow. */

static inline int stencil_index (scalar s, IJK i)
{
  int len = 1, index = 0;
  for (int d = 0; d < dimension; d++) {
    if (i.i[d] < 0 || i.i[d] >= _STENCIL) // stencil overflow
      return -1;
    index += len*i.i[d], len *= _STENCIL;
  }
  return index;
}

/**
This is the reverse of the function above. */

static inline IJK stencil_ijk (int index)
{
  IJK i;
  int len = _STENCIL_SIZE;
  for (int d = dimension - 1; d >= 0; d--) {
    len /= _STENCIL;
    i.i[d] = index/len;
    index -= len*i.i[d];
    i.i[d] -= _STENCIL/2;
  }
  return i;
}

/**
This displays a (1D,2D,3D) stencil index. */

static void write_stencil_index (IJK i, int shift)
{
  sysfprintf (qstderr(), "[%d", i.i[0] - shift);
  for (int d = 1; d < dimension; d++)
    sysfprintf (qstderr(), ",%d", i.i[d] - shift);
  sysfprintf (qstderr(), "]");
}

/**
This function is called whenever a field index `i` is called
(i.e. typically when `s[i,j,k]` code is executed). It checks for
stencil overflows and updates the `read` stencil array. */

int _stencil_access (scalar s, IJK i, const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    return 0;
  int index = stencil_index (s, i);
  if (index < 0) {
    sysfprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
		file, line, s.name);
    write_stencil_index (i, _STENCIL/2);
    sysfprintf (qstderr(), "\n");
    fflush (qstderr());
    abort();
  }
  s.read[index]++;
  return index;
}

/**
The `foreach_stencil()` loop definition replaces the foreach(),
foreach_face() and foreach_vertex() loops. It implements the
read/write access detection used for the definition of stencils.

Before loop execution the `read` and `write` arrays are
initialised. Non-trivial values are set in the `write` array. */

@def foreach_stencil() {
  for (scalar _s in baseblock) {
    for (int _i = 0; _i < _STENCIL_SIZE; _i++) {
      _s.read[_i] = 0;
      _s.write[_i] = 1.7759437274 + _s.i + _i;
    }
  }
  _foreach_data.face = 0;
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  // fixme: this assumes that Point is always something like:
  // struct { int i,.., k, level, n, ...; }
  // with i,..,k the dimension indices, level the level and
  // n the number of grid points
  for (int _i = 0; _i < dimension; _i++)
    ((int *)&point)[_i] = _STENCIL/2;
  if (sizeof(Point) >= (dimension + 2)*sizeof(int))
    ((int *)&point)[dimension + 1] = 1;
  POINT_VARIABLES
@

/**
After the stencil loop execution the `write` array is checked for any
modification. */
    
@def end_foreach_stencil()
  for (scalar _s in baseblock) {
    for (int _i = 0; _i < _STENCIL_SIZE; _i++)
      _s.write[_i] = (_s.write[_i] != 1.7759437274 + _s.i + _i);
  }
  end_stencil (&_foreach_data);
}
@

@define foreach_face_stencil foreach_stencil
@define end_foreach_face_stencil() end_foreach_stencil()

@define foreach_vertex_stencil foreach_stencil
@define end_foreach_vertex_stencil() end_foreach_stencil()

@define is_stencil_face_x() ((_foreach_data.face |= (1 << 0)))
@define is_stencil_face_y() ((_foreach_data.face |= (1 << 1)))
@define is_stencil_face_z() ((_foreach_data.face |= (1 << 2)))

/**
This is the function called in case of missing reductions. */

void reduction_warning (const char * fname, int line, const char * var)
{
  fprintf (stderr,
  "%s:%d: warning: variable '%s' is modified by this foreach loop:\n"
  "%s:%d: warning: use a loop-local variable, a reduction operation\n"
  "%s:%d: warning: or a serial loop to get rid of this warning\n",
	   fname, line, var, fname, line, fname, line);
}

/**
## Automatic boundary conditions

Boundary conditions need to be applied if `s` is dirty, or if any of
the field `d` it depends on is dirty. */

static inline bool scalar_is_dirty (scalar s)
{
  if (s.dirty)
    return true;
  scalar * depends = s.depends;
  for (scalar d in depends)
    if (d.dirty)
      return true;
  return false;
}

/**
Does the boundary conditions on `a` depend on those on `b`? */

static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = a.depends;
  for (scalar s in depends)
    if (s.i == b.i)
      return true;
  return false;
}

/**
There are two types of boundary conditions: "full" boundary
conditions, done by `boundary_internal()` and "flux" boundary
conditions (i.e. normal components on faces only) done by
`boundary_face()`. */

void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face)  (vectorl);

/**
This function is called after the stencil access detection, just
before the (real) foreach loop is executed. This is where we use the
stencil access pattern to see whether boundary conditions need to be
applied. */

void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  loop->vertex = !strcmp (loop->each, "foreach_vertex");
  
  /**
  We check the accesses for each field... */
  
  for (scalar s in baseblock) {
    bool write = false, read = false;
    int max = 0; // maximum stencil width/height/depth

    /**
    ... and for each stencil position. */
    
    for (int n = 0; n < _STENCIL_SIZE; n++)
      if (s.write[n] || s.read[n]) {
	IJK i = stencil_ijk (n);

	/**
	This stencil position has been write-accessed. Only the
	central index ([0,0,0]) can be write-accessed. */
	
	if (s.write[n]) {
	  for (int d = 0; d < dimension; d++)
	    if (i.i[d] != 0) {
	      fprintf (stderr,
		       "%s:%d: error: illegal write within this loop: %s",
		       loop->fname, loop->line, s.name);
	      write_stencil_index (i, 0);
	      fprintf (stderr, "\n");
	      fflush (stderr);
	      abort();
	    }
	  write = true;
	}

	/**
	This stencil position has been read-accessed, we record the
	maximum "width" of the stencil. */
	
	if (s.read[n]) {
	  read = true; // fixme: what if write == true and read == 1 ??
	  int d = 0;
	  foreach_dimension() {
	    if ((!s.face || s.v.x.i != s.i) && abs(i.i[d]) > max)
	      max = abs(i.i[d]);
	    d++;
	  }
	}
      }
    
    /**
    We now know whether this field has been read- and/or write-accessed. */
    
#ifdef foreach_layer
    if (_layer == 0 || s.block == 1)
#endif
    {

      /**
      If the field is read and dirty, we need to check if boundary
      conditions need to be applied. */
      
      if (read && scalar_is_dirty (s)) {

	/**
	If this is a face field, we check whether "full" BCs need to
	be applied, or whether "face" BCs are sufficient. */
	
	if (s.face) {
	  if (max > 0) // face, stencil wider than 0
	    listc = list_append (listc, s);
	  else if (!write) { // face, flux only
	    scalar sn = s.v.x.i >= 0 ? s.v.x : s;
	    foreach_dimension()
	      if (s.v.x.i == s.i) {

		/* fixme: imposing BCs on fluxes should be done by
		   boundary_face() .*/
		
		if (sn.boundary[left] || sn.boundary[right])
		  listc = list_append (listc, s);
		else if (s.dirty != 2) {
		  listf.x = list_append (listf.x, s);
		  flux = true;
		}
	      }
	  }
	}

	/**
	For dirty, centered fields BCs need to be applied if the
	stencil is wider than zero. */
	
	else if (max > 0)
	  listc = list_append (listc, s);
      }

      /**
      Write accesses need to be consistent with the declared field
      type (i.e. face or vertex). */
      
      if (write) {
	if (dimension > 1 && !loop->vertex && loop->first) {
	  bool vertex = true;
	  foreach_dimension()
	    if (s.d.x != -1)
	      vertex = false;
	  if (vertex)
	    fprintf (stderr,
		     "%s:%d: warning: vertex scalar '%s' should be assigned with"
		     " a foreach_vertex() loop\n",
		     loop->fname, loop->line, s.name);
	}
	if (s.face) {
	  if (loop->face == 0 && loop->first)
	    fprintf (stderr,
		     "%s:%d: warning: face vector '%s' should be assigned with"
		     " a foreach_face() loop\n",
		     loop->fname, loop->line, s.name);
	}
	else if (loop->face) {
	  if (s.v.x.i < 0) { // scalar
	    int d = 1, i = 0;
	    foreach_dimension() {
	      if (loop->face == d) {
		s.face = 2, s.v.x.i = s.i;
		s.boundary[left] = s.boundary[right] = NULL;
#if PRINTBOUNDARY
		fprintf (stderr,
			 "%s:%d: turned %s into a face vector %c-component\n",
			 loop->fname, loop->line, s.name, 'x' + i);
#endif
	      }
	      d *= 2, i++;
	    }
	    if (!s.face && loop->first)
	      fprintf (stderr,
		       "%s:%d: warning: scalar '%s' should be assigned with "
		       "a foreach_face(x|y|z) loop\n",
		       loop->fname, loop->line, s.name);
	  }
	  else { // vector
	    char * name = NULL;
	    if (s.name) {
	      name = strdup (s.name);
	      char * s = name + strlen(name) - 1;
	      while (s != name && *s != '.') s--;
	      if (s != name) *s = '\0';
	    }
	    init_face_vector (s.v, name);
#if PRINTBOUNDARY
	    fprintf (stderr, "%s:%d: turned %s into a face vector\n",
		     loop->fname, loop->line, name);
#endif
	    free (name);
	  }
	}
	else if (loop->vertex) {
	  bool vertex = true;
	  foreach_dimension()
	    if (s.d.x != -1)
	      vertex = false;
	  if (!vertex) {
	    char * name = NULL;
	    if (s.name) name = strdup (s.name);
	    init_vertex_scalar (s, name);
	    foreach_dimension()
	      s.v.x.i = -1;
#if PRINTBOUNDARY
	    fprintf (stderr, "%s:%d: turned %s into a vertex scalar\n",
		     loop->fname, loop->line, name);
#endif
	    free (name);
	  }
	}

	/**
	If the field is write-accessed, we add it to the 'dirty'
	list. */
	
	dirty = list_append (dirty, s);
	for (scalar d in baseblock)
	  if (scalar_depends_from (d, s))
	    dirty = list_append (dirty, d);
      }
    }
  }

  /**
  We apply face (flux) boundary conditions. */
  
  if (flux) {
#if PRINTBOUNDARY
    int i = 0;
    foreach_dimension() {
      if (listf.x) {
	fprintf (stderr, "%s:%d: flux %c:", loop->fname, loop->line, 'x' + i);
	for (scalar s in listf.x)
	  fprintf (stderr, " %d:%s", s.i, s.name);
	fputc ('\n', stderr);
      }
      i++;
    }
#endif
    boundary_face (listf);
    foreach_dimension()
      free (listf.x);
  }
  
  /**
  We apply "full" boundary conditions. */

  if (listc) {
#if PRINTBOUNDARY
    fprintf (stderr, "%s:%d: listc:", loop->fname, loop->line);
    for (scalar s in listc)
      fprintf (stderr, " %d:%s", s.i, s.name);
    fputc ('\n', stderr);
#endif
    boundary_internal (listc, loop->fname, loop->line);
    free (listc);
  }

  /**
  We update the dirty status of fields which will be write-accessed by
  the foreach loop. */
  
  if (dirty) {
#if PRINTBOUNDARY
    fprintf (stderr, "%s:%d: dirty:", loop->fname, loop->line);
    for (scalar s in dirty)
      fprintf (stderr, " %d:%s", s.i, s.name);
    fputc ('\n', stderr);
#endif
    for (scalar s in dirty)
      s.dirty = true;
    free (dirty);
  }
}

/**
## See also

* [Stencil test case](/src/test/stencils.c)
*/
