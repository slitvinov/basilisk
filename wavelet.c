#ifndef foreach_fine_to_coarse
#  error "the grid needs to implement foreach_fine_to_coarse()"
#endif

void restriction (Data * m, int n, var a)
{
  foreach_fine_to_coarse (m, n) {
    coarse(a,0,0) = (fine(a,0,0) + fine(a,1,0) + fine(a,0,1) + fine(a,1,1))/4.;
  } end_foreach_fine_to_coarse();
}

void wavelet (Data * m, int n, var a, var w)
{
  foreach_fine_to_coarse (m, n) {
    /* difference between fine value and bilinearly-interpolated coarse value */
    fine(w,0,0) = fine(a,0,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,-1)) + coarse(a,-1,-1))/16.;
    fine(w,0,1) = fine(a,0,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,-1,0) + coarse(a,0,+1)) + coarse(a,-1,+1))/16.;
    fine(w,1,0) = fine(a,1,0) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,0,-1) + coarse(a,+1,0)) + coarse(a,+1,-1))/16.;
    fine(w,1,1) = fine(a,1,1) - 
      (9.*coarse(a,0,0) + 3.*(coarse(a,+1,0) + coarse(a,0,+1)) + coarse(a,+1,+1))/16.;
  } end_foreach_fine_to_coarse ();
}
