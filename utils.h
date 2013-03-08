#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>

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
  fflush (fp);
}

void output_matrix (scalar f, int n, FILE * fp, bool linear)
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
      float y = delta*j - 0.5 + delta/2., v;
      if (linear)
	v = interpolate (f, x, y);
      else {
	Point point = locate (x, y);
	VARIABLES;
	v = val (f, 0, 0);
      }
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
}

void output_cells (FILE * fp)
{
  foreach() {
    delta /= 2.;
    fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	     x - delta, y - delta,
	     x - delta, y + delta,
	     x + delta, y + delta,
	     x + delta, y - delta,
	     x - delta, y - delta);
  }
}

typedef struct {
  clock_t c;
  struct timeval tv;
} timer_t;

timer_t timer_start (void)
{
  timer_t t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
  return t;
}

void timer_print (timer_t t, int i, int tnc)
{
  clock_t end = clock ();
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  double cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  double real = (tvend.tv_sec - t.tv.tv_sec) + (tvend.tv_usec - t.tv.tv_usec)/1e6;
  if (tnc < 0) {
    tnc = 0;
    foreach(reduction(+:tnc)) tnc++;
    tnc *= i;
  }
  fprintf (stderr, "# " GRIDNAME ", %d steps, %g CPU, %.4g real, %.3g points.steps/s\n",
	   i, cpu, real, tnc/real);
}

typedef struct {
  double avg, rms, max, area;
} norm;

norm normf (scalar f)
{
  norm n = { 0., 0., 0., 0. };
  foreach(reduction(max:n.max) reduction(+:n.avg) reduction(+:n.rms) reduction(+:n.area)) {
    double v = fabs(f[]);
    if (v > n.max) n.max = v;
    double a = sq(delta);
    n.avg  += a*v;
    n.rms  += a*v*v;
    n.area += a;
  }
  n.avg /= n.area;
  n.rms = sqrt(n.rms/n.area);
  return n;
}

typedef struct {
  double min, max, sum;
} stats;

stats statsf (scalar f)
{
  stats s = { 1e100, -1e100, 0. };
  foreach(reduction(+:s.sum) reduction(max:s.max) reduction(min:s.min)) {
    s.sum += f[]*delta*delta;
    if (f[] > s.max) s.max = f[];
    if (f[] < s.min) s.min = f[];
  }
  return s;
}
