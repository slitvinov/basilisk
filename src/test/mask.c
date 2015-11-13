/**
Checks that "half-mapped" fine/coarse cells do not cause trouble with
restriction of face values. 

This only works when two layers of ghost cells are used on the
boundaries. */

#define BGHOSTS 2
#define TRASH 1
#include "grid/quadtree.h"

int main() {
  face vector v[];
  mask (y >  L0/4. ? top : none);
  foreach_face()
    v.x[] = 1.;
  boundary ((scalar *){v});
  restriction ((scalar *){v});  
}
