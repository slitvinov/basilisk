#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#define GHOSTS  1 // number of ghost layers
#define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

#define pi 3.14159265358979
#undef HUGE
#define HUGE 1e30
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }

#ifdef _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#else
# define OMP(x)
#endif
// fixme: _OMPSTART and _OMPEND are only used for working around the
// lack of min|max reduction operations in OpenMP < 3.1
#define OMP_PARALLEL()     OMP(omp parallel) { _OMPSTART
#define OMP_END_PARALLEL() _OMPEND }
#define _OMPSTART
#define _OMPEND

#define NOT_UNUSED(x) (x = x)

#define VARIABLES

void * grid = NULL;       // the grid

typedef int scalar;

typedef struct {
  scalar x, y;
} vector;

typedef struct {
  vector x, y;
} tensor;

#define val(a,k,l)     data(k,l)[a]
#define fine(a,k,l)    _fine(a,k,l)
#define coarse(a,k,l)  _coarse(a,k,l)
#define neighbor(k,l)  _neighbor(k,l)

#define scalars(...) (scalar []){__VA_ARGS__,-1}
#define vectors(...) (vector []){__VA_ARGS__,{-1,-1}}
#define none         (scalar []){-1}

int vectors_len (vector * list)
{
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

int scalars_len (scalar * list)
{
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * scalars_append (scalar * list, scalar a)
{
  int ns = scalars_len (list);
  scalar * list1 = malloc ((ns + 2)*sizeof (scalar));
  ns = 0;
  for (scalar s in list) {
    assert (s != a); // a is already in the list
    list1[ns++] = s;
  }
  list1[ns++] = a;
  list1[ns] = -1;
  return list1;
}

// basic methods

scalar (* new_scalar) (scalar);
vector (* new_vector) (vector);
tensor (* new_tensor) (tensor);

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
typedef double (* BoundaryFunc) (Point, scalar);
typedef void   (* RefineFunc)   (Point, scalar);

BoundaryFunc * _boundary[nboundary]; // boundary conditions for each scalar
RefineFunc   * _refine;              // refinement function for each scalar

void free_solver()
{
  for (int b = 0; b < nboundary; b++)
    free (_boundary[b]);
  free (_refine);
}
