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
// CFL number
double CFL = 0.5;

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
  if (point.level < 0)
    return nodata;
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

void output_field (scalar * list, int n, FILE * fp, bool linear)
{
  fprintf (fp, "# 1:x 2:y");
  int i = 3;
  for (scalar s in list)
    fprintf (fp, " %d:%d", i++, s);
  fputc('\n', fp);
  double delta = L0/n;
  for (int i = 0; i < n; i++) {
    double x = delta*i + X0 + delta/2.;
    for (int j = 0; j < n; j++) {
      double y = delta*j + Y0 + delta/2.;
      fprintf (fp, "%g %g", x, y);
      if (linear) {
	for (scalar s in list)
	  fprintf (fp, " %g", interpolate (s, x, y));
      }
      else {
	Point point = locate (x, y);
	for (scalar s in list)
	  fprintf (fp, " %g", point.level >= 0 ? val(s,0,0) : nodata);
      }
      fputc ('\n', fp);
    }
    fputc ('\n', fp);
  }
  fflush (fp);
}

void output_matrix (scalar f, int n, FILE * fp, bool linear)
{
  float fn = n;
  float delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float y = delta*j + X0 + delta/2.;
    fwrite (&y, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float x = delta*i + X0 + delta/2.;
    fwrite (&x, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float y = delta*j + Y0 + delta/2., v;
      if (linear)
	v = interpolate (f, x, y);
      else {
	Point point = locate (x, y);
	assert (point.level >= 0);
	v = val(f,0,0);
      }
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
}

typedef struct {
  clock_t c;
  struct timeval tv;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
  return t;
}

void timer_print (timer t, int i, int tnc)
{
  clock_t end = clock ();
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  double cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  double real = ((tvend.tv_sec - t.tv.tv_sec) + 
		 (tvend.tv_usec - t.tv.tv_usec)/1e6);
  if (tnc < 0) {
    tnc = 0;
    foreach(reduction(+:tnc)) tnc++;
    tnc *= i;
  }
  printf ("# " GRIDNAME 
	  ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %ld var\n",
	  i, cpu, real, tnc/real, datasize/sizeof(double));
}

typedef struct {
  double avg, rms, max, area;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., area = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:area)) {
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

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}

double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}

double zero (double s0, double s1, double s2)
{
  return 0.;
}

void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  trash ((scalar *) g);
  foreach() {
    scalar s; vector v;
    for (s,v in f,g)
      if (s.gradient)
	foreach_dimension()
	  v.x[] = s.gradient (s[-1,0], s[], s[1,0])/delta;
      else // centered
	foreach_dimension()
	  v.x[] = (s[1,0] - s[-1,0])/(2.*delta);
  }
  boundary ((scalar *) g);
}
