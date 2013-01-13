#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#if 0

typedef struct {
  double u[3][3], v[3][3], h[3][3], b[3][3];
  double un, vn, hn;
} Data;
// #define INDEX(a,k,l) m[i][j].a[k-i+1][l-j+1]
#define INDEX(a,k,l) m[i+k][j+l].a[0][0]

#else

typedef struct {
  double u, v, h, b;
  double un, vn, hn;
} Data;
#define INDEX(a,k,l) m[i+k][j+l].a

#endif

#define u(k,l) INDEX(u,k,l)
#define v(k,l) INDEX(v,k,l)
#define h(k,l) INDEX(h,k,l)
#define b(k,l) INDEX(b,k,l)
#define un(k,l) m[i+k][j+l].un
#define vn(k,l) m[i+k][j+l].vn
#define hn(k,l) m[i+k][j+l].hn

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

void update_u (Data ** m, int n)
{
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      u(0,0) = un(0,0);
      v(0,0) = vn(0,0);
    }
}

void update_h (Data ** m, int n)
{
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      h(0,0) = hn(0,0);
}

void boundary_h (Data ** m, int n)
{
  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    h(0,-1) = h(0,0);
    j = n;
    h(0,1) = h(0,0);
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    h(-1,0) = h(0,0);
    i = n;
    h(1,0) = h(0,0);
  }
}

void boundary_u (Data ** m, int n)
{
  int i, j;
  for (i = 1; i <= n; i++) {
    /* solid walls */
    j = 1; v(0,0) = 0.;
    j = n; v(0,1) = 0.;
  }
  for (j = 1; j <= n; j++) {
    /* solid walls */
    i = 1; u(0,0) = 0.;
    i = n; u(1,0) = 0.;
  }
}

void tracer_advection (Data ** m, int n, double dt)
{
  double h = L0/n;
  dt /= 2.*h;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      hn(0,0) = h(0,0) + dt*((h(0,0) + h(-1,0))*u(0,0) - 
			     (h(0,0) + h(1,0))*u(1,0) +
			     (h(0,0) + h(0,-1))*v(0,0) - 
			     (h(0,0) + h(0,1))*v(0,1));
}

void tracer_advection_upwind (Data ** m, int n, double dt)
{
  double h = L0/n;
  dt /= h;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      hn(0,0) = h(0,0) + dt*((u(0,0) < 0. ? h(0,0) : h(-1,0))*u(0,0) - 
			     (u(1,0) > 0. ? h(0,0) : h(1,0))*u(1,0) +
			     (v(0,0) < 0. ? h(0,0) : h(0,-1))*v(0,0) - 
			     (v(0,1) > 0. ? h(0,0) : h(0,1))*v(0,1));
}

static double vorticity (Data ** m, int i, int j, int n)
{
  return (v(0,0) - v(-1,0) + u(0,-1) - u(0,0))/(L0/n);
}

static double KE (Data ** m, int i, int j)
{
#if 1
  double uc = u(0,0) + u(1,0);
  double vc = v(0,0) + v(0,1);
  return (uc*uc + vc*vc)/8.;
#else
  double uc = u(0,0)*u(0,0) + u(1,0)*u(1,0);
  double vc = v(0,0)*v(0,0) + v(0,1)*v(0,1);
  return (uc + vc)/4.;
#endif
}

void momentum (Data ** m, int n, double dt)
{
  double dtg = dt*G/(L0/n);
  double dtf = dt/4.;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      double vort = vorticity(m,i,j,n);
      double g = h(0,0) + b(0,0) + KE(m,i,j);
      double vortu = (vort + vorticity(m,i,j+1,n))/2.;
      un(0,0) = u(0,0)
	- dtg*(g - h(-1,0) - b(-1,0) - KE(m,i-1,j))
	+ dtf*(vortu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1));
      double vortv = (vort + vorticity (m,i+1,j,n))/2.;
      vn(0,0) = v(0,0)
	- dtg*(g - h(0,-1) - b(0,-1) - KE(m,i,j-1))
	- dtf*(vortv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1));
    }
}

#include "init.h"

void output_field (Data ** m, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      double x = XC, y = YC;
      fprintf (fp, "%g %g %g\n", x, y, h(0,0) + b(0,0));
    }
    fprintf (fp, "\n");
  }
}

double timestep (Data ** m, int n)
{
  double dx = L0/n;
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  dx *= dx;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) {
      if (h(0,0) > 0.) {
	double dt = dx/(G*h(0,0));
	if (dt < dtmax) dtmax = dt;
      }
      if (u(0,0) != 0.) {
	double dt = dx/(u(0,0)*u(0,0));
	if (dt < dtmax) dtmax = dt;
      }
      if (v(0,0) != 0.) {
	double dt = dx/(v(0,0)*v(0,0));
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
  
  Data ** m = matrix_new (n + 2, n + 2, sizeof (Data));

  initial_conditions (m, n);
  boundary_h (m, n);
  boundary_u (m, n);

  clock_t start, end;
  start = clock ();
  do {
    #include "output.h"
    double dt = timestep (m, n);
    tracer_advection (m, n, dt);
    update_h (m, n);
    boundary_h (m, n);
    momentum (m, n, dt);
    update_u (m, n);
    boundary_u (m, n);
    t += dt; i++;
  } while (t < TMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
