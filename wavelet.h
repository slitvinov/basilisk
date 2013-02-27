#ifndef foreach_fine_to_coarse
#  error "the grid needs to implement foreach_fine_to_coarse()"
#endif

void restriction (void * grid, scalar start, scalar end)
{
  foreach_fine_to_coarse (grid)
    for (scalar v = start; v <= end; v++)
      val(v,0,0) = (fine(v,0,0) + fine(v,1,0) + fine(v,0,1) + fine(v,1,1))/4.;
}

void restriction_u_v (void * grid, scalar u, scalar v)
{
  foreach_fine_to_coarse (grid) {
    u[] = (fine(u,0,0) + fine(u,0,1))/2.;
    v[] = (fine(v,0,0) + fine(v,1,0))/2.;
  }
}

void wavelet (void * grid, scalar v, scalar w)
{
  foreach_fine_to_coarse (grid) {
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
