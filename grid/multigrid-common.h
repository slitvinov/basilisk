#include "cartesian-common.h"

void restriction (scalar * list)
{
  foreach_fine_to_coarse()
    for (scalar s in list)
      s[] = (fine(s,0,0) + fine(s,1,0) + fine(s,0,1) + fine(s,1,1))/4.;
}

void wavelet (scalar v, scalar w)
{
  foreach_fine_to_coarse() {
    /* difference between fine value and bilinearly-interpolated coarse value */
    fine(w,0,0) = fine(v,0,0) - 
      (9.*v[] + 3.*(v[-1,0] + v[0,-1]) + v[-1,-1])/16.;
    fine(w,0,1) = fine(v,0,1) - 
      (9.*v[] + 3.*(v[-1,0] + v[0,+1]) + v[-1,+1])/16.;
    fine(w,1,0) = fine(v,1,0) - 
      (9.*v[] + 3.*(v[+1,0] + v[0,-1]) + v[+1,-1])/16.;
    fine(w,1,1) = fine(v,1,1) - 
      (9.*v[] + 3.*(v[+1,0] + v[0,+1]) + v[+1,+1])/16.;
  }
}
