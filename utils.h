#include <stdio.h>
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
#define VARIABLES					       \
  double delta = DELTA;          /* cell size */	       \
  double x  = XC(I), y  = YC(J); /* cell center */	       \
  double xu = XU(I), yu = YU(J); /* staggered u-coordinates */ \
  double xv = XV(I), yv = YV(J); /* staggered v-coordinates */ \
  /* we need this to avoid compiler warnings */	               \
  NOT_UNUSED(x);  NOT_UNUSED(y);        \
  NOT_UNUSED(xu); NOT_UNUSED(yu);	\
  NOT_UNUSED(xv); NOT_UNUSED(yv);	\
  NOT_UNUSED(delta);

void symmetry (scalar v)
{
  foreach_boundary (right)  v[+1,0] = v[];
  foreach_boundary (left)   v[-1,0] = v[];
  foreach_boundary (top)    v[0,+1] = v[];
  foreach_boundary (bottom) v[0,-1] = v[];
}

void uv_symmetry (scalar u, scalar v)
{
  foreach_boundary (right) {
    u[1,0] = 0.;
    v[1,0] = v[];
  }
  foreach_boundary (left) {
    u[-1,0] = - u[1,0];
    u[] = 0.;
    v[-1,0] = v[];
  }
  foreach_boundary (top) {
    v[0,1] = 0.;
    u[0,1] = u[];
  }
  foreach_boundary (bottom) {
    v[0,-1] = - v[0,1];
    v[] = 0.;
    u[0,-1] = u[];
  }
}

void runge_kutta (int stages,
		  double t, double dt,
		  int nv, scalar f[nv], scalar df[stages][nv], 
		  void (* advance) (double t, scalar f[nv], scalar df[nv]),
		  void (* update)  (double t, scalar f[nv]))
{
  switch (stages) {
  case 1:
    (* advance) (t, f, df[0]);
    foreach()
      for (int v = 0; v < nv; v++)
	f[v][] += df[0][v][]*dt;
    (* update) (t + dt, f);
    break;

  case 2:
    (* advance) (t, f, df[0]);
    foreach()
      for (int v = 0; v < nv; v++)
	df[0][v][] = f[v][] + df[0][v][]*dt/2.;
    (* update) (t + dt/2., df[0]);

    (* advance) (t + dt/2., df[0], df[1]);
    foreach()
      for (int v = 0; v < nv; v++)
	f[v][] += df[1][v][]*dt;
    (* update) (t + dt, f);
    break;

  default:
    /* not implemented yet */
    assert(false);
  }
}

double change (scalar v, scalar vn)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    double dv = fabs (v[] - vn[]);
    if (dv > max)
      max = dv;
    vn[] = v[];
  }
  return max;
}

double interpolate (scalar v, double xp, double yp)
{
  Point point = locate (xp, yp);
  VARIABLES;
  x = (xp - x)/delta;
  y = (yp - y)/delta;
  assert (x >= -0.5 && x <= 0.5);
  assert (y >= -0.5 && y <= 0.5);
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);
  /* bilinear interpolation */
  return (val(v,0,0)*(1. - x)*(1. - y) + 
	  val(v,i,0)*x*(1. - y) + 
	  val(v,0,j)*(1. - x)*y + 
	  val(v,i,j)*x*y);
}

void output_field (scalar f, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  double delta = 1./n;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double x = delta*i - 0.5 + delta/2., y = delta*j - 0.5 + delta/2.;
      fprintf (fp, "%g %g %g\n", x, y, interpolate (f, x, y));
    }
    fputc ('\n', fp);
  }
}

void output_matrix (scalar f, int n, FILE * fp)
{
  float fn = n;
  float delta = 1./fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float y = delta*j - 0.5 + delta/2.;
    fwrite (&y, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float x = delta*i - 0.5 + delta/2.;
    fwrite (&x, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float y = delta*j - 0.5 + delta/2.;
      float v = interpolate (f, x, y);
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
}
