#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// Default parameters, do not change them!! edit parameters.h instead
// number of grid points
int N = 64;
// maximum timestep
double DT = 1e10;
// size of the box
double L0 = 1.;
// acceleration of gravity
double G = 1.;
// Coriolis parameter
double F0 = 1.;
// end time
double TMAX = 1e10;
// CFL number
double CFL = 0.5;

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

void symmetry_conditions (double ** t, int n)
{
  for (int i = 1; i <= n; i++) {
    t[i][0] = t[i][1];
    t[i][n+1] = t[i][n];
    t[0][i] = t[1][i];
    t[n+1][i] = t[n][i];
  }      
}

void solid_walls_conditions (double ** u, double ** v, int n)
{
  for (int i = 1; i <= n; i++)
    u[1][i] = u[n+1][i] = v[i][1] = v[i][n+1] = 0.;
}

void tracer_advection (double ** t, double ** u, double ** v, int n, double dt,
		       double ** tn)
{
  double h = L0/n;
  dt /= 2.*h;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      tn[i][j] = t[i][j] + dt*((t[i][j] + t[i-1][j])*u[i][j] - 
			       (t[i][j] + t[i+1][j])*u[i+1][j] +
			       (t[i][j] + t[i][j-1])*v[i][j] - 
			       (t[i][j] + t[i][j+1])*v[i][j+1]);
  symmetry_conditions (tn, n);
}

void tracer_advection_upwind (double ** t, double ** u, double ** v, int n, double dt,
			      double ** tn)
{
  double h = L0/n;
  dt /= h;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      tn[i][j] = t[i][j] + dt*((u[i][j] < 0. ? t[i][j] : t[i-1][j])*u[i][j] - 
			       (u[i+1][j] > 0. ? t[i][j] : t[i+1][j])*u[i+1][j] +
			       (v[i][j] < 0. ? t[i][j] : t[i][j-1])*v[i][j] - 
			       (v[i][j+1] > 0. ? t[i][j] : t[i][j+1])*v[i][j+1]);
  symmetry_conditions (tn, n);
}

static double vorticity (double ** u, double ** v, int i, int j, int n)
{
  return (v[i][j] - v[i-1][j] + u[i][j-1] - u[i][j])/(L0/n);
}

static double KE (double ** u, double ** v, int i, int j)
{
#if 1
  double uc = u[i][j] + u[i+1][j];
  double vc = v[i][j] + v[i][j+1];
  return (uc*uc + vc*vc)/8.;
#else
  double uc = u[i][j]*u[i][j] + u[i+1][j]*u[i+1][j];
  double vc = v[i][j]*v[i][j] + v[i][j+1]*v[i][j+1];
  return (uc + vc)/4.;
#endif
}

void momentum (double ** h, double ** b,
	       double ** u, double ** v, int n, double dt,
	       double ** un, double ** vn)
{
  double dtg = dt*G/(L0/n);
  double dtf = dt/4.;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      double vort = vorticity(u,v,i,j,n);
      double g = h[i][j] + b[i][j] + KE(u,v,i,j);
      double vortu = (vort + vorticity(u,v,i,j+1,n))/2.;
      un[i][j] = u[i][j]
	- dtg*(g - h[i-1][j] - b[i-1][j] - KE(u,v,i-1,j))
	+ dtf*(vortu + F0)*(v[i][j] + v[i][j+1] + v[i-1][j] + v[i-1][j+1]);
      double vortv = (vort + vorticity (u, v, i+1, j, n))/2.;
      vn[i][j] = v[i][j]
	- dtg*(g - h[i][j-1] - b[i][j-1] - KE(u,v,i,j-1))
	- dtf*(vortv + F0)*(u[i][j] + u[i+1][j] + u[i][j-1] + u[i+1][j-1]);
    }
  solid_walls_conditions (un, vn, n);
}

#include "init.h"

void output_field (double ** h, double ** b, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      double x = XC, y = YC;
      fprintf (fp, "%g %g %g\n", x, y, h[i][j] + b[i][j]);
    }
    fprintf (fp, "\n");
  }
}

#define swap(a,b) {double ** tmp = a; a = b; b = tmp;}

double timestep (double ** u, double ** v, double ** h, int n)
{
  double dx = L0/n;
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  dx *= dx;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      if (h[i][j] > 0.) {
	double dt = dx/(G*h[i][j]);
	if (dt < dtmax) dtmax = dt;
      }
      if (u[i][j] != 0.) {
	double dt = dx/(u[i][j]*u[i][j]);
	if (dt < dtmax) dtmax = dt;
      }
      if (v[i][j] != 0.) {
	double dt = dx/(v[i][j]*v[i][j]);
	if (dt < dtmax) dtmax = dt;
      }
    }
  return sqrt (dtmax)*CFL;
}

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = N;
  double ** u = matrix_new (n + 2, n + 2, sizeof (double));
  double ** v = matrix_new (n + 2, n + 2, sizeof (double));
  double ** h = matrix_new (n + 2, n + 2, sizeof (double));
  double ** b = matrix_new (n + 2, n + 2, sizeof (double));
  
  double ** un = matrix_new (n + 2, n + 2, sizeof (double));
  double ** vn = matrix_new (n + 2, n + 2, sizeof (double));
  double ** hn = matrix_new (n + 2, n + 2, sizeof (double));

  initial_conditions (u, v, h, b, n);

  clock_t start, end;
  start = clock ();
  do {
    #include "output.h"
    double dt = timestep (u, v, h, n);
    tracer_advection (h, u, v, n, dt, hn);
    swap (h, hn);
    momentum (h, b, u, v, n, dt, un, vn);
    swap (u, un);
    swap (v, vn);
    t += dt; i++;
  } while (t < TMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "%d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
