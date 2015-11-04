/**
Checks that "half-mapped" fine/coarse cells do not cause trouble with
restriction of face values. */

#define TRASH 1
#include "grid/quadtree.h"

int main() {
  face vector v[];
  mask (y >  L0/4. ? top : none);
  foreach_face()
    v.x[] = 0.;
  boundary ((scalar *){v});
  restriction ((scalar *){v});  
}
