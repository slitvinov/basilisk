#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#define GHOSTS  1 // number of ghost layers
#define TRASH   1 // whether to 'trash' uninitialised data 
                  // (useful for debugging)

#define pi 3.14159265358979
#define undefined 1e100
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

scalar new_scalar (scalar s);
vector new_vector (vector v);
tensor new_tensor (tensor t);

#define scalars(...) (scalar []){__VA_ARGS__,-1}
#define vectors(...) (vector []){__VA_ARGS__,{-1,-1}}
#define val(a,k,l) data(k,l)[a]
#define fine(a,k,l)    _fine(a,k,l)
#define coarse(a,k,l)  _coarse(a,k,l)
#define neighbor(k,l)  _neighbor(k,l)

// events 

typedef struct _Event Event;
typedef int (* Expr) (int *, double *);

struct _Event {
  int last, nexpr;
  int (* action) (int, double);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  double t;
  int i, a;
};

double tnext = undefined; // time of next event
void init_events (void);

// boundary conditions for each direction/variable

enum { right, left, top, bottom, nboundary };
// ghost cell coordinates for each direction
int _ig[nboundary] = {1,-1,0,0}, 
    _jg[nboundary] = {0,0,1,-1};

void ** _boundary[nboundary];

void free_boundaries ()
{
  for (int b = 0; b < nboundary; b++)
    free (_boundary[b]);
}
