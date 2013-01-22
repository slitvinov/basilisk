#ifndef foreach_fine_to_coarse
#  error "the grid needs to implement foreach_fine_to_coarse()"
#endif

void restriction (void * grid, var a)
{
  foreach_fine_to_coarse (grid)
    val(a,0,0) = (fine(a,0,0) + fine(a,1,0) + fine(a,0,1) + fine(a,1,1))/4.;
}

void wavelet (void * grid, var a, var w)
{
  foreach_fine_to_coarse (grid) {
    /* difference between fine value and bilinearly-interpolated coarse value */
    fine(w,0,0) = fine(a,0,0) - 
      (9.*val(a,0,0) + 3.*(val(a,-1,0) + val(a,0,-1)) + val(a,-1,-1))/16.;
    fine(w,0,1) = fine(a,0,1) - 
      (9.*val(a,0,0) + 3.*(val(a,-1,0) + val(a,0,+1)) + val(a,-1,+1))/16.;
    fine(w,1,0) = fine(a,1,0) - 
      (9.*val(a,0,0) + 3.*(val(a,+1,0) + val(a,0,-1)) + val(a,+1,-1))/16.;
    fine(w,1,1) = fine(a,1,1) - 
      (9.*val(a,0,0) + 3.*(val(a,+1,0) + val(a,0,+1)) + val(a,+1,+1))/16.;
  }
}
