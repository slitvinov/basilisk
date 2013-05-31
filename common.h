#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#define GHOSTS  1 // number of ghost layers
#define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

#define pi 3.14159265358979
#undef HUGE
#define HUGE 1e30
#define nodata DBL_MAX
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }

#ifdef _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
# define pid() omp_get_thread_num()
#else
# define OMP(x)
# define pid() 0
#endif
// fixme: _OMPSTART and _OMPEND are only used for working around the
// lack of min|max reduction operations in OpenMP < 3.1
#define OMP_PARALLEL()     OMP(omp parallel) { _OMPSTART
#define OMP_END_PARALLEL() _OMPEND }
#define _OMPSTART
#define _OMPEND

#define NOT_UNUSED(x) (x = x)

#define VARIABLES

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

#define norm(v) (sqrt(sq(val(v.x,0,0)) + sq(val(v.y,0,0))))

#define val(a,k,l)     data(k,l)[a]
#define fine(a,k,l)    _fine(a,k,l)
#define coarse(a,k,l)  _coarse(a,k,l)
#define neighbor(k,l)  _neighbor(k,l)

// lists

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = realloc (list, sizeof (vector)*(len + 2));
  list[len] = v;
  list[len + 1] = (vector){-1,-1};
  return list;
}

int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
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

scalar * all = NULL; // all the fields

// basic methods

scalar (* init_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);

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
  double t;
  int i, a;
};

double tnext = INFINITY; // time of next event
void init_events (void);

// boundary conditions for each direction/variable

enum { right, left, top, bottom, nboundary };
// ghost cell coordinates for each direction
int _ig[nboundary] = {1,-1,0,0}, 
    _jg[nboundary] = {0,0,1,-1};

typedef struct _Point Point;

// methods for each scalar

typedef struct {
  double (* boundary[nboundary]) (Point, scalar);
  void   (* refine)              (Point, scalar);
  void   (* coarsen)             (Point, scalar);
  double (* gradient)            (double, double, double);
} Methods;

Methods * _method;

void free_solver()
{
  free (_method); _method = NULL;
  free (all); all = NULL;
  grid = NULL;
}
