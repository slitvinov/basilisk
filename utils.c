#include <math.h>

// Default parameters, do not change them!! edit parameters.h instead
// number of grid points
int N = 64;
// maximum timestep
double DT = 1e10;
// coordinates of the center of the box
double X0 = -0.5, Y0 = -0.5;
// size of the box
double L0 = 1.;
// end time
double TMAX = 1e10;
// end iterations
int IMAX = 1 << 30;
// CFL number
double CFL = 0.5;

void symmetry (Data * m, int n, var v)
{
  foreach_boundary (m, n, right)  { stencil(v,+1,0) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, left)   { stencil(v,-1,0) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, top)    { stencil(v,0,+1) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, bottom) { stencil(v,0,-1) = val(v); } end_foreach_boundary();  
}

void uv_symmetry (Data * m, int n, var u, var v)
{
  foreach_boundary (m, n, right) {
    stencil(u,+1,0) = 0.;
    stencil(v,+1,0) = val(v);
  } end_foreach_boundary();
  foreach_boundary (m, n, left) {
    stencil(u,-1,0) = 0.;
    stencil(v,-1,0) = val(v);
  } end_foreach_boundary();
  foreach_boundary (m, n, top) {
    stencil(v,0,+1) = 0.;
    stencil(u,0,+1) = val(u);
  } end_foreach_boundary();
  foreach_boundary (m, n, bottom) {
    stencil(v,0,-1) = 0.;
    stencil(u,0,-1) = val(u);
  } end_foreach_boundary();
}

void * matrix_new (int n, int p, int size)
{
  int i;
  void ** m;
  char * a;
  
  m = malloc (n*sizeof (void **));
  a = malloc (n*p*size);
  for (i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return (void *) m;
}

void matrix_free (void * m)
{
  free (((void **) m)[0]);
  free (m);
}

