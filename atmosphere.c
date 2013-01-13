#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define M(i,j) m[(i)*(n + 2) + (j)]
#if 0
  #define NS 3
  #define INDEX(a,k,l) M(i,j).a[k+NS/2][l+NS/2]
//  #define INDEX(a,k,l) M(i+k,j+l).a[NS/2][NS/2]
#else
  #define NS 1
  #define INDEX(a,k,l) M(i+k,j+l).a[NS/2][NS/2]
#endif

typedef struct {
  double u[NS][NS], v[NS][NS], h[NS][NS], b[NS][NS];
  double un, vn, hn;
} Data;

#define u(k,l) INDEX(u,k,l)
#define v(k,l) INDEX(v,k,l)
#define h(k,l) INDEX(h,k,l)
#define b(k,l) INDEX(b,k,l)
#define un(k,l) M(i+k,j+l).un
#define vn(k,l) M(i+k,j+l).vn
#define hn(k,l) M(i+k,j+l).hn

#define foreach(m) for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++)

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

#include "boundary.h"

void tracer_advection (Data * m, int n, double dt)
{
  double h = L0/n;
  dt /= 2.*h;
  foreach (m)
    hn(0,0) = h(0,0) + dt*((h(0,0) + h(-1,0))*u(0,0) - 
			   (h(0,0) + h(1,0))*u(1,0) +
			   (h(0,0) + h(0,-1))*v(0,0) - 
			   (h(0,0) + h(0,1))*v(0,1));
}

void tracer_advection_upwind (Data * m, int n, double dt)
{
  double h = L0/n;
  dt /= h;
  foreach (m)
    hn(0,0) = h(0,0) + dt*((u(0,0) < 0. ? h(0,0) : h(-1,0))*u(0,0) - 
			   (u(1,0) > 0. ? h(0,0) : h(1,0))*u(1,0) +
			   (v(0,0) < 0. ? h(0,0) : h(0,-1))*v(0,0) - 
			   (v(0,1) > 0. ? h(0,0) : h(0,1))*v(0,1));
}

static double vorticity (Data * m, int i, int j, int n)
{
  return (v(0,0) - v(-1,0) + u(0,-1) - u(0,0))/(L0/n);
}

static double KE (Data * m, int i, int j, int n)
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

void momentum (Data * m, int n, double dt)
{
  double dtg = dt*G/(L0/n);
  double dtf = dt/4.;
  foreach (m) {
    double vort = vorticity(m,i,j,n);
    double g = h(0,0) + b(0,0) + KE(m,i,j,n);
    double vortu = (vort + vorticity(m,i,j+1,n))/2.;
    un(0,0) = u(0,0)
      - dtg*(g - h(-1,0) - b(-1,0) - KE(m,i-1,j,n))
      + dtf*(vortu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1));
    double vortv = (vort + vorticity (m,i+1,j,n))/2.;
    vn(0,0) = v(0,0)
      - dtg*(g - h(0,-1) - b(0,-1) - KE(m,i,j-1,n))
      - dtf*(vortv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1));
  }
}

#include "init.h"

void output_field (Data * m, int n, FILE * fp)
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

double timestep (Data * m, int n)
{
  double dx = L0/n;
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  dx *= dx;
  foreach (m) {
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
  
  Data * m = malloc ((n + 2)*(n + 2)*sizeof (Data));

  initial_conditions (m, n);
  boundary_b (m, n);
  boundary_h (m, n);
  boundary_u (m, n);

  clock_t start, end;
  start = clock ();
  do {
    #include "output.h"
    double dt = timestep (m, n);
    tracer_advection (m, n, dt);
    boundary_h (m, n);
    momentum (m, n, dt);
    boundary_u (m, n);
    t += dt; i++;
  } while (t < TMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
