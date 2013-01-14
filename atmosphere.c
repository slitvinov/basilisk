#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// Default parameters, do not change them!! edit parameters.h instead
// number of grid points
int N = 64;
// maximum timestep
double DT = 1e10;
// coordinates of the center of the box
double X0 = -0.5, Y0 = -0.5;
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

#define DX (L0/n)
#define XC(i) ((i + 0.5)*DX + X0)
#define YC(j) ((j + 0.5)*DX + X0)
#define XU(i) ((i)*DX + X0)
#define YU(j) YC(j)
#define XV(i) XC(i)
#define YV(j) ((j)*DX + X0)

#include "grid.h"

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

void tracer_advection (Data * m, int n, double dt)
{
  dt /= 2.*DX;
  foreach (m)
    hn(0,0) = h(0,0) + dt*((h(0,0) + h(-1,0))*u(0,0) - 
			   (h(0,0) + h(1,0))*u(1,0) +
			   (h(0,0) + h(0,-1))*v(0,0) - 
			   (h(0,0) + h(0,1))*v(0,1));
}

void tracer_advection_upwind (Data * m, int n, double dt)
{
  dt /= DX;
  foreach (m)
    hn(0,0) = h(0,0) + dt*((u(0,0) < 0. ? h(0,0) : h(-1,0))*u(0,0) - 
			   (u(1,0) > 0. ? h(0,0) : h(1,0))*u(1,0) +
			   (v(0,0) < 0. ? h(0,0) : h(0,-1))*v(0,0) - 
			   (v(0,1) > 0. ? h(0,0) : h(0,1))*v(0,1));
}

void momentum (Data * m, int n, double dt)
{
  double dtg = dt*G/DX;
  double dtf = dt/4.;
  foreach (m) {
    double g = h(0,0) + b(0,0) + ke(0,0);
    double psiu = (psi(0,0) + psi(0,1))/2.;
    un(0,0) = u(0,0)
      - dtg*(g - h(-1,0) - b(-1,0) - ke(-1,0))
      + dtf*(psiu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1));
    double psiv = (psi(0,0) + psi(1,0))/2.;
    vn(0,0) = v(0,0)
      - dtg*(g - h(0,-1) - b(0,-1) - ke(0,-1))
      - dtf*(psiv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1));
  }
}

void ke_psi (Data * m, int n)
{
  foreach (m) {
#if 1
    double uc = u(0,0) + u(1,0);
    double vc = v(0,0) + v(0,1);
    ke(0,0) = (uc*uc + vc*vc)/8.;
#else
    double uc = u(0,0)*u(0,0) + u(1,0)*u(1,0);
    double vc = v(0,0)*v(0,0) + v(0,1)*v(0,1);
    ke(0,0) = (uc + vc)/4.;
#endif
    psi(0,0) = (v(0,0) - v(-1,0) + u(0,-1) - u(0,0))/DX;
  }    
}

double timestep (Data * m, int n)
{
  double dx = DX;
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

#include "init.h"

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = N;

  Data * m = init_grid (n);
  initial_conditions (m, n);
  boundary_b (m, n);
  boundary_h (m, n);
  boundary_u (m, n);
  ke_psi (m, n);
  boundary_ke_psi (m, n);

  clock_t start, end;
  start = clock ();
  do {
    double dt = timestep (m, n);
    #include "output.h"
    tracer_advection (m, n, dt);
    boundary_h (m, n);
    momentum (m, n, dt);
    boundary_u (m, n);
    ke_psi (m, n);
    boundary_ke_psi (m, n);
    t += dt; i++;
  } while (t < TMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));
}
