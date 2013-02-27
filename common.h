#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#define GHOSTS 1  // number of ghost layers
#define TRASH  1  // whether to 'trash' uninitialised data (useful for debugging)

#define NOT_UNUSED(x) (x = x)

#define VARIABLES

enum { right, left, top, bottom };

typedef int scalar;
typedef struct {
  scalar x, y;
} vector;

#define val(a,k,l) data(k,l)[a]

#define pi 3.14159265358979
#define undefined 1e100
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)

typedef struct _Event Event;
typedef int (* Expr) (int *, double *);

struct _Event {
  int last, nexpr;
  void (* action) (int, double);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  double t;
  int i, a;
};

void * grid = NULL;       // the grid
double tnext = undefined; // time of next event
