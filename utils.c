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

#define DX    (L0*delta)
#define XC(i) ((i + 0.5)*DX + X0)
#define YC(j) ((j + 0.5)*DX + X0)
#define XU(i) ((i)*DX + X0)
#define YU(j) YC(j)
#define XV(i) XC(i)
#define YV(j) ((j)*DX + X0)

#undef VARIABLES
#define VARIABLES \
  double x  = XC(I), y  = YC(J); /* cell center */	\
  double xu = XU(I), yu = YU(J); /* staggered u-coordinates */ \
  double xv = XV(I), yv = YV(J); /* staggered v-coordinates */ \
  /* we need this to avoid compiler warnings */	\
  NOT_UNUSED(x);  NOT_UNUSED(y);		        \
  NOT_UNUSED(xu); NOT_UNUSED(yu);                 \
  NOT_UNUSED(xv); NOT_UNUSED(yv);

void symmetry (void * grid, var v)
{
  foreach_boundary (grid, right)  { val(v,+1,0) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (grid, left)   { val(v,-1,0) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (grid, top)    { val(v,0,+1) = val(v,0,0); } end_foreach_boundary();
  foreach_boundary (grid, bottom) { val(v,0,-1) = val(v,0,0); } end_foreach_boundary();  
}

void uv_symmetry (void * grid, var u, var v)
{
  foreach_boundary (grid, right) {
    val(u,+1,0) = 0.;
    val(v,+1,0) = val(v,0,0);
  } end_foreach_boundary();
  foreach_boundary (grid, left) {
    val(u,-1,0) = 0.;
    val(v,-1,0) = val(v,0,0);
  } end_foreach_boundary();
  foreach_boundary (grid, top) {
    val(v,0,+1) = 0.;
    val(u,0,+1) = val(u,0,0);
  } end_foreach_boundary();
  foreach_boundary (grid, bottom) {
    val(v,0,-1) = 0.;
    val(u,0,-1) = val(u,0,0);
  } end_foreach_boundary();
}

void runge_kutta (int stages, double dt,
		  void * grid, 
		  int nv, var f[nv], var df[stages][nv], 
		  void (* advance) (void * grid, var f[nv], var df[nv]),
		  void (* update)  (void * grid, var f[nv]))
{
  switch (stages) {
  case 1:
    (* advance) (grid, f, df[0]);
    foreach (grid)
      for (int v = 0; v < nv; v++)
	val(f[v],0,0) = val(f[v],0,0) + val(df[0][v],0,0)*dt;
    end_foreach();
    (* update) (grid, f);
    break;

  case 2:
    (* advance) (grid, f, df[0]);
    foreach (grid)
      for (int v = 0; v < nv; v++)
	val(df[0][v],0,0) = val(f[v],0,0) + val(df[0][v],0,0)*dt/2.;
    end_foreach();
    (* update) (grid, df[0]);

    (* advance) (grid, df[0], df[1]);
    foreach (grid)
      for (int v = 0; v < nv; v++)
	val(f[v],0,0) += val(df[1][v],0,0)*dt;
    end_foreach();
    (* update) (grid, f);
    break;

  default:
    /* not implemented yet */
    assert(false);
  }
}
