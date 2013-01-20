#define GRID "Cartesian grid"

#include "utils.h"
#include <stdlib.h>

#define I (i - 1)
#define J (j - 1)

#define data(k,l) _m[(i + k)*(_n + 2) + (j + l)]

#define foreach(m,n) {				\
  Data * _m = m; int _n = n;			\
  for (int i = 1; i <= _n; i++)			\
    for (int j = 1; j <= _n; j++) { \
      VARIABLES
#define end_foreach() }}

#define foreach_boundary(m,n,d) {		 \
  Data * _m = m; int _n = n;			 \
  for (int _k = 1; _k <= _n; _k++) {		 \
    int i = d > left ? _k : d == right ? _n : 1; \
    int j = d < top  ? _k : d == top   ? _n : 1; \
    VARIABLES
#define end_foreach_boundary() }}

void * init_grid (int n)
{
  return calloc ((n + 2)*(n + 2), sizeof (Data));
}

#define free_grid free
