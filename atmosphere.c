#include <time.h>

struct _Data {
  double u, v, h, b, ke, psi;
  double un, vn, hn;
};

#define u(k,l)    data(k,l).u
#define v(k,l)    data(k,l).v
#define h(k,l)    data(k,l).h
#define b(k,l)    data(k,l).b
#define ke(k,l)   data(k,l).ke
#define psi(k,l)  data(k,l).psi
#define un(k,l)   data(k,l).un
#define vn(k,l)   data(k,l).vn
#define hn(k,l)   data(k,l).hn

#include "utils.h"
#include "grid.c"
#include "utils.c"

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;

void tracer_advection (Data * m, int n, double dt)
{
  dt /= 2.*(L0/n); /* fixme */
  foreach (m, n)
    hn(0,0) = h(0,0) + dt*((h(0,0) + h(-1,0))*u(0,0) - 
			   (h(0,0) + h(1,0))*u(1,0) +
			   (h(0,0) + h(0,-1))*v(0,0) - 
			   (h(0,0) + h(0,1))*v(0,1));
  end_foreach();
}

void tracer_advection_upwind (Data * m, int n, double dt)
{
  dt /= (L0/n); /* fixme */
  foreach (m, n)
    hn(0,0) = h(0,0) + dt*((u(0,0) < 0. ? h(0,0) : h(-1,0))*u(0,0) - 
			   (u(1,0) > 0. ? h(0,0) : h(1,0))*u(1,0) +
			   (v(0,0) < 0. ? h(0,0) : h(0,-1))*v(0,0) - 
			   (v(0,1) > 0. ? h(0,0) : h(0,1))*v(0,1));
  end_foreach();
}

double timestep (Data * m, int n)
{
  double dx = (L0/n); /* fixme */
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  dx *= dx;
  foreach (m, n) {
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
  } end_foreach();
  return sqrt (dtmax)*CFL;
}

void momentum (Data * m, int n, double dt)
{
  double dtg = dt*G/(L0/n); /* fixme */
  double dtf = dt/4.;
  foreach (m, n) {
    double g = h(0,0) + b(0,0) + ke(0,0);
    double psiu = (psi(0,0) + psi(0,1))/2.;
    un(0,0) = u(0,0)
      - dtg*(g - h(-1,0) - b(-1,0) - ke(-1,0))
      + dtf*(psiu + F0)*(v(0,0) + v(0,1) + v(-1,0) + v(-1,1));
    double psiv = (psi(0,0) + psi(1,0))/2.;
    vn(0,0) = v(0,0)
      - dtg*(g - h(0,-1) - b(0,-1) - ke(0,-1))
      - dtf*(psiv + F0)*(u(0,0) + u(1,0) + u(0,-1) + u(1,-1));
  } end_foreach();
}

void ke_psi (Data * m, int n)
{
  foreach (m, n) {
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
  } end_foreach();
}

#include "init.h"

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = N;

  void * m = init_grid (n);
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
    foreach (m, n) { h(0,0) = hn(0,0); } end_foreach();
    boundary_h (m, n);
    momentum (m, n, dt);
    foreach (m, n) {
      u(0,0) = un(0,0);
      v(0,0) = vn(0,0);
    } end_foreach();
    boundary_u (m, n);
    ke_psi (m, n);
    boundary_ke_psi (m, n);
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRID ", %d timesteps, %g CPU, %d points.steps/s\n",
	   i, cpu, (int) (n*n*i/cpu));

  free_grid (m);
}
