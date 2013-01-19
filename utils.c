#include <math.h>
#include <assert.h>

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

void symmetry (void * m, int n, var v)
{
  foreach_boundary (m, n, right)  { val(v,+1,0) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (m, n, left)   { val(v,-1,0) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (m, n, top)    { val(v,0,+1) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (m, n, bottom) { val(v,0,-1) = val(v,0,0); } end_foreach_boundary();  
}

void uv_symmetry (void * m, int n, var u, var v)
{
  foreach_boundary (m, n, right) {
    val(u,+1,0) = 0.;
    val(v,+1,0) = val(v,0,0);
  } end_foreach_boundary();
  foreach_boundary (m, n, left) {
    val(u,-1,0) = 0.;
    val(v,-1,0) = val(v,0,0);
  } end_foreach_boundary();
  foreach_boundary (m, n, top) {
    val(v,0,+1) = 0.;
    val(u,0,+1) = val(u,0,0);
  } end_foreach_boundary();
  foreach_boundary (m, n, bottom) {
    val(v,0,-1) = 0.;
    val(u,0,-1) = val(u,0,0);
  } end_foreach_boundary();
}

void runge_kutta (int stages, double dt,
		  void * m, int n, 
		  int nv, var f[nv], var df[stages][nv], 
		  void (* advance) (void * m, int n, var f[nv], var df[nv]),
		  void (* update)  (void * m, int n, var f[nv]))
{
  switch (stages) {
  case 1:
    (* advance) (m, n, f, df[0]);
    foreach (m, n)
      for (int v = 0; v < nv; v++)
	val(f[v],0,0) = val(f[v],0,0) + val(df[0][v],0,0)*dt;
    end_foreach();
    (* update) (m, n, f);
    break;

  case 2:
    (* advance) (m, n, f, df[0]);
    foreach (m, n)
      for (int v = 0; v < nv; v++)
	val(df[0][v],0,0) = val(f[v],0,0) + val(df[0][v],0,0)*dt/2.;
    end_foreach();
    (* update) (m, n, df[0]);

    (* advance) (m, n, df[0], df[1]);
    foreach (m, n)
      for (int v = 0; v < nv; v++)
	val(f[v],0,0) += val(df[1][v],0,0)*dt;
    end_foreach();
    (* update) (m, n, f);
    break;

  default:
    /* not implemented yet */
    assert(false);
  }
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
