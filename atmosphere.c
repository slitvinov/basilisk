#include "utils.h"

struct _Data {
  int flags; /* fixme: this is grid-specific */
  double u, v, h, b, ke, psi;
  double un, vn, hn;
};

#define u(k,l)    celln(k,l).u
#define v(k,l)    celln(k,l).v
#define h(k,l)    celln(k,l).h
#define b(k,l)    celln(k,l).b
#define ke(k,l)   celln(k,l).ke
#define psi(k,l)  celln(k,l).psi
#define un(k,l)   celln(k,l).un
#define vn(k,l)   celln(k,l).vn
#define hn(k,l)   celln(k,l).hn

#include "grid.h"
#include "utils.c"

#include <time.h>

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
}
