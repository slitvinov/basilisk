@include <stdlib.h>
@include <stdio.h>
@include <stddef.h>
@include <stdbool.h>
@include <string.h>
@include <float.h>
@include <limits.h>
@include <assert.h>
@include <math.h>

@define pi 3.14159265358979
@undef HUGE
@define HUGE ((double)1e30)
@define nodata DBL_MAX

@define max(a,b) ((a) > (b) ? (a) : (b))
@define min(a,b) ((a) < (b) ? (a) : (b))
@define sq(x) ((x)*(x))
@define sign(x) ((x) > 0 ? 1 : -1)
@define noise() (1. - 2.*rand()/(double)RAND_MAX)
@define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
@define swap(type,a,b) { type tmp = a; a = b; b = tmp; }

@define GHOSTS  1 // number of ghost layers
@define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

@if _OPENMP
@ include <omp.h>
@ define OMP(x) Pragma(#x)
@ define pid() omp_get_thread_num()
@else
@ define OMP(x)
@ define pid() 0
@endif
// fixme: _OMPSTART and _OMPEND are only used for working around the
// lack of min|max reduction operations in OpenMP < 3.1
@define OMP_PARALLEL()     OMP(omp parallel) { _OMPSTART
@define OMP_END_PARALLEL() _OMPEND }
@define _OMPSTART
@define _OMPEND

@define NOT_UNUSED(x) (x = x)

@define VARIABLES
@define val(a,k,l)     data(k,l)[a]
@define fine(a,k,l)    _fine(a,k,l)
@define coarse(a,k,l)  _coarse(a,k,l)
@define allocated(k,l) _allocated(k,l)
@define neighbor(k,l)  _neighbor(k,l)

// the grid
void * grid = NULL;
// coordinates of the lower-left corner of the box
double X0 = 0., Y0 = 0.;
// size of the box
double L0 = 1.;

typedef int scalar;

typedef struct {
  scalar x, y;
} vector;

typedef struct {
  vector x, y;
} tensor;

#define norm(v) (sqrt(sq(v.x[]) + sq(v.y[])))

// lists

int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t == s)
      return list;
  int len = list_len (list);
  list = realloc (list, sizeof (scalar)*(len + 2));
  list[len] = s;
  list[len + 1] = -1;
  return list;
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    for (scalar s1 in l)
      if (s1 == s)
	return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    for (scalar s in l)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  for (scalar s in l2)
    l3 = list_append (l3, s);
  return l3;
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  for (vector w in list)
    if (w.x == v.x && w.y == v.y)
      return list;
  int len = vectors_len (list);
  list = realloc (list, sizeof (vector)*(len + 2));
  list[len] = v;
  list[len + 1] = (vector){-1,-1};
  return list;
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    for (vector v in l)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (*s >= 0) {
    vector v;
    v.x = *s++;
    assert (*s >= 0);
    v.y = *s++;
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  for (tensor t in list) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = realloc (list, sizeof (tensor)*(len + 2));
  list[len] = t;
  list[len + 1] = (tensor){{-1,-1},{-1,-1}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x >= 0) {
    tensor t;
    t.x = *v++;
    assert (v->x >= 0);
    t.y = *v++;
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL; // all the fields

// basic methods

scalar (* init_scalar)           (scalar, const char *);
vector (* init_vector)           (vector, const char *);
tensor (* init_tensor)           (tensor, const char *);
vector (* init_staggered_vector) (vector, const char *);

// events 

typedef struct _Event Event;
typedef int (* Expr) (int *, double *);

struct _Event {
  int last, nexpr;
  int (* action) (int, double);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
};

double tnext = HUGE; // time of next event
void init_events (void);

// boundary conditions for each direction/variable

enum { right, left, top, bottom, nboundary };
// ghost cell coordinates for each direction
int _ig[nboundary] = {1,-1,0,0}, 
    _jg[nboundary] = {0,0,1,-1};

@define dirichlet(x)            (2.*(x) - val(_s,0,0))
@define dirichlet_homogeneous() (- val(_s,0,0))
@define neumann(x)              (Delta*(x) + val(_s,0,0))
@define neumann_homogeneous()   (val(_s,0,0))

typedef struct _Point Point;

// methods for each scalar

typedef struct {
  double (* boundary[nboundary])             (Point, scalar);
  double (* boundary_homogeneous[nboundary]) (Point, scalar);
  void   (* refine)                          (Point, scalar);
  void   (* coarsen)                         (Point, scalar);
  double (* gradient)                        (double, double, double);
  struct { int x, y; } d;        // staggering
  vector v;
  bool   staggered;
} Methods;

Methods * _method;

void free_solver()
{
  free (_method); _method = NULL;
  free (all); all = NULL;
  grid = NULL;
}
