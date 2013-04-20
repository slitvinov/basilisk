#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

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

#define GHOSTS  1 // number of ghost layers
#define TRASH   1 // whether to 'trash' uninitialised data (useful for debugging)
#define TWO_ONE 1 // enforce 2:1 refinement ratio

#define NOT_UNUSED(x) (x = x)

#define VARIABLES

enum { right, left, top, bottom, nboundary };
// ghost cell coordinates for each direction
int _ig[nboundary] = {1,-1,0,0}, 
    _jg[nboundary] = {0,0,1,-1};

typedef int scalar;

typedef struct {
  scalar x, y;
} vector;

typedef struct {
  vector x, y;
} tensor;

#define val(a,k,l) data(k,l)[a]

#define pi 3.14159265358979
#define undefined 1e100
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }

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
void * grid = NULL;       // the grid
typedef void (* Boundary) (int l);

Boundary * _boundary[nboundary]; // boundary conditions for each direction/variable

void init_boundaries (int nvar)
{
  for (int b = 0; b < nboundary; b++)
    _boundary[b] = calloc (nvar, sizeof (Boundary));
}

void free_boundaries ()
{
  for (int b = 0; b < nboundary; b++)
    free (_boundary[b]);
}
