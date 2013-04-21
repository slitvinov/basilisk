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
#define XC(i) ((i + 0.5)*DX + X0*L0)
#define YC(j) ((j + 0.5)*DX + X0*L0)
#define XU(i) ((i)*DX + X0*L0)
#define YU(j) YC(j)
#define XV(i) XC(i)
#define YV(j) ((j)*DX + X0*L0)

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
  printf ("# " GRIDNAME ", %d steps, %g CPU, %.4g real, %.3g points.steps/s\n",
	  i, cpu, real, tnc/real);
}

typedef struct {
  double avg, rms, max, area;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., area = 0.;
  foreach(reduction(max:max) reduction(+:avg) reduction(+:rms) reduction(+:area)) {
    double v = fabs(f[]);
    if (v > max) max = v;
    double a = sq(delta);
    avg  += a*v;
    rms  += a*v*v;
    area += a;
  }
  norm n;
  n.avg = avg/area;
  n.rms = sqrt(rms/area);
  n.max = max;
  n.area = area;
  return n;
}

typedef struct {
  double min, max, sum;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0.;
  foreach(reduction(+:sum) reduction(max:max) reduction(min:min)) {
    sum += f[]*delta*delta;
    if (f[] > max) max = f[];
    if (f[] < min) min = f[];
  }
  stats s;
  s.min = min, s.max = max, s.sum = sum;
  return s;
}

static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

static double minmod1   (double r) { return generic_limiter (r, 1.);  }
static double superbee1 (double r) { return generic_limiter (r, 2.);  }
static double sweby1    (double r) { return generic_limiter (r, 1.5); }

void centered (const scalar f, vector g)
{
  foreach()
    foreach_dimension()
      g.x[] = (f[1,0] - f[-1,0])/2./delta;
}

void minmod (const scalar f, vector g)
{
  foreach()
    foreach_dimension()
      g.x[] = minmod1 ((f[1,0] - f[])/(f[] - f[-1,0]))*(f[] - f[-1,0])/delta;
}

double theta = 1.3;

static double minmod2 (double a, double b, double c)
{
  if (a > 0. && b > 0. && c > 0.)
    return min(min(a,b),c);
  if (a < 0. && b < 0. && c < 0.)
    return max(max(a,b),c);
  return 0.;
}

void generalized_minmod (scalar * f, vector * g)
{
  /* see (A.6) in 
   *    Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
   *    Saint-Venant system. Mathematical Modelling and Numerical
   *    Analysis, 36(3), 397-425.*/
  foreach() {
    scalar s; vector v;
    for (s,v in f,g)
      foreach_dimension()
	v.x[] = minmod2 (theta*(s[] - s[-1,0]), 
			 (s[1,0] - s[-1,0])/2., 
			 theta*(s[1,0] - s[]))/delta;
  }
}

void superbee (scalar * f, vector * g)
{
  foreach() {
    scalar s; vector v;
    for (s,v in f,g)
      foreach_dimension()
	v.x[] = superbee1 ((s[1,0] - s[])/(s[] - s[-1,0]))*(s[] - s[-1,0])
	/delta;
  }
}

void sweby (scalar * f, vector * g)
{
  foreach() {
    scalar s; vector v;
    for (s,v in f,g)
      foreach_dimension()
	v.x[] = sweby1 ((s[1,0] - s[])/(s[] - s[-1,0]))*(s[] - s[-1,0])/delta;
  }
}

void zero (scalar * f, vector * g)
{
  foreach()
    for (vector v in g)
      foreach_dimension()
	v.x[] = 0.;
}
