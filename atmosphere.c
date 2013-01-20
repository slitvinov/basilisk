#include <time.h>

struct _Data {
  double u, v, h, b, ke, psi;
  double un, vn, hn;
  //  double un1, vn1, hn1;
};

#define h(k,l)    data(k,l).h
#define u(k,l)    data(k,l).u
#define v(k,l)    data(k,l).v

#define b(k,l)    data(k,l).b
#define ke(k,l)   data(k,l).ke
#define psi(k,l)  data(k,l).psi
#define un(k,l)   data(k,l).un
#define vn(k,l)   data(k,l).vn
#define hn(k,l)   data(k,l).hn

#include "grid.c"
#include "utils.c"

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;

void advection_centered (void * m, int n, var f, var u, var v, var df)
{
  double dx = L0/n; // fixme
  foreach (m, n)
    val(df,0,0) = ((val(f,0,0) + val(f,-1,0))*val(u,0,0) - 
		   (val(f,0,0) + val(f,1,0))*val(u,1,0) +
		   (val(f,0,0) + val(f,0,-1))*val(v,0,0) - 
		   (val(f,0,0) + val(f,0,1))*val(v,0,1))/(2.*dx);
  end_foreach();
}

void advection_upwind (void * m, int n, var f, var u, var v, var df)
{
  double dx = L0/n; // fixme
  foreach (m, n)
    val(df,0,0) = ((val(u,0,0) < 0. ? val(f,0,0) : val(f,-1,0))*val(u,0,0) - 
		   (val(u,1,0) > 0. ? val(f,0,0) : val(f,1,0))*val(u,1,0) +
		   (val(v,0,0) < 0. ? val(f,0,0) : val(f,0,-1))*val(v,0,0) - 
		   (val(v,0,1) > 0. ? val(f,0,0) : val(f,0,1))*val(v,0,1))/dx;
  end_foreach();
}

double timestep (void * m, int n)
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

void momentum (void * m, int n, var u, var v, var h, var du, var dv)
{
  double dx = L0/n; // fixme
  foreach (m, n) {
    double g = val(h,0,0) + b(0,0) + ke(0,0);
    double psiu = (psi(0,0) + psi(0,1))/2.;
    val(du,0,0) = 
      - G*(g - val(h,-1,0) - b(-1,0) - ke(-1,0))/dx
      + (psiu + F0)*(val(v,0,0) + val(v,0,1) + val(v,-1,0) + val(v,-1,1))/4.;
    double psiv = (psi(0,0) + psi(1,0))/2.;
    val(dv,0,0) = 
      - G*(g - val(h,0,-1) - b(0,-1) - ke(0,-1))/dx
      - (psiv + F0)*(val(u,0,0) + val(u,1,0) + val(u,0,-1) + val(u,1,-1))/4.;
  } end_foreach();
}

void ke_psi (void * m, int n, var u, var v)
{
  foreach (m, n) {
#if 1
    double uc = val(u,0,0) + val(u,1,0);
    double vc = val(v,0,0) + val(v,0,1);
    ke(0,0) = (uc*uc + vc*vc)/8.;
#else
    double uc = val(u,0,0)*val(u,0,0) + val(u,1,0)*val(u,1,0);
    double vc = val(v,0,0)*val(v,0,0) + val(v,0,1)*val(v,0,1);
    ke(0,0) = (uc + vc)/4.;
#endif
    psi(0,0) = (val(v,0,0) - val(v,-1,0) + val(u,0,-1) - val(u,0,0))/DX;
  } end_foreach();
}

void advance (void * m, int n, var * f, var * df)
{
  var u = f[0], v = f[1], h = f[2];
  var du = df[0], dv = df[1], dh = df[2];

  advection_centered (m, n, h, u, v, dh);
  momentum (m, n, u, v, h, du, dv);
}

#include "init.h"

void update (void * m, int n, var * f)
{
  var u = f[0], v = f[1], h = f[2];
  boundary_h (m, n, h);
  boundary_u (m, n, u, v);
  ke_psi (m, n, u, v);
  boundary_ke_psi (m, n);
}

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0, n = N;

  void * m = init_grid (n);
  initial_conditions (m, n);
  boundary_b (m, n);
  boundary_h (m, n, var(h));
  boundary_u (m, n, var(u), var(v));
  ke_psi (m, n, var(u), var(v));
  boundary_ke_psi (m, n);

  clock_t start, end;
  start = clock ();
  do {
    double dt = timestep (m, n);
    #include "output.h"
#if 1
    advection_centered (m, n, var(h), var(u), var(v), var(hn));
    foreach (m, n) { h(0,0) += hn(0,0)*dt; } end_foreach();
    boundary_h (m, n, var(h));
    momentum (m, n, var(u), var(v), var(h), var(un), var(vn));
    foreach (m, n) {
      u(0,0) += un(0,0)*dt;
      v(0,0) += vn(0,0)*dt;
    } end_foreach();
    boundary_u (m, n, var(u), var(v));
    ke_psi (m, n, var(u), var(v));
    boundary_ke_psi (m, n);
#else /* unstable! */
    var f[3] =  { var(u), var(v), var(h) };
    var df[2][3] = {{ var(un),  var(vn),  var(hn) },
		    { var(un1), var(vn1), var(hn1) }};
    runge_kutta (2, dt, m, n, 3, f, df, advance, update);
#endif
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRID ", %d timesteps, %g CPU, %g points.steps/s\n",
	   i, cpu, (n*n*(double)i/cpu));

  free_grid (m);
}
