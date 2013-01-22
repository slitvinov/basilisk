#include <time.h>

struct _Data {
  double u, v, h, b, ke, psi;
  double un, vn, hn;
  double h1, e;
  //  double un1, vn1, hn1;
};

#include GRID
#include "utils.h"

var u = var(u), v = var(v), h = var(h), b = var(b), ke = var(ke), psi = var(psi),
  un = var(un), vn = var(vn), hn = var(hn);

// Default parameters, do not change them!! edit parameters.h instead
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;

void advection_centered (void * grid, var f, var u, var v, var df)
{
  foreach (grid)
    val(df,0,0) = ((val(f,0,0) + val(f,-1,0))*val(u,0,0) - 
		   (val(f,0,0) + val(f,1,0))*val(u,1,0) +
		   (val(f,0,0) + val(f,0,-1))*val(v,0,0) - 
		   (val(f,0,0) + val(f,0,1))*val(v,0,1))/(2.*L0*delta);
}

void advection_upwind (void * grid, var f, var u, var v, var df)
{
  foreach (grid)
    val(df,0,0) = ((val(u,0,0) < 0. ? val(f,0,0) : val(f,-1,0))*val(u,0,0) - 
		   (val(u,1,0) > 0. ? val(f,0,0) : val(f,1,0))*val(u,1,0) +
		   (val(v,0,0) < 0. ? val(f,0,0) : val(f,0,-1))*val(v,0,0) - 
		   (val(v,0,1) > 0. ? val(f,0,0) : val(f,0,1))*val(v,0,1))/(L0*delta);
}

double timestep (void * grid)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach (grid) {
    double dx = L0*delta;
    dx *= dx;
    if (val(h,0,0) > 0.) {
      double dt = dx/(G*val(h,0,0));
      if (dt < dtmax) dtmax = dt;
    }
    if (val(u,0,0) != 0.) {
      double dt = dx/(val(u,0,0)*val(u,0,0));
      if (dt < dtmax) dtmax = dt;
    }
    if (val(v,0,0) != 0.) {
      double dt = dx/(val(v,0,0)*val(v,0,0));
      if (dt < dtmax) dtmax = dt;
    }
  }
  return sqrt (dtmax)*CFL;
}

void momentum (void * grid, var u, var v, var h, var du, var dv)
{
  foreach (grid) {
    double g = val(h,0,0) + val(b,0,0) + val(ke,0,0);
    double psiu = (val(psi,0,0) + val(psi,0,1))/2.;
    double dx = L0*delta;
    val(du,0,0) = 
      - G*(g - val(h,-1,0) - val(b,-1,0) - val(ke,-1,0))/dx
      + (psiu + F0)*(val(v,0,0) + val(v,0,1) + val(v,-1,0) + val(v,-1,1))/4.;
    double psiv = (val(psi,0,0) + val(psi,1,0))/2.;
    val(dv,0,0) = 
      - G*(g - val(h,0,-1) - val(b,0,-1) - val(ke,0,-1))/dx
      - (psiv + F0)*(val(u,0,0) + val(u,1,0) + val(u,0,-1) + val(u,1,-1))/4.;
  }
}

void ke_psi (void * grid, var u, var v)
{
  foreach (grid) {
#if 1
    double uc = val(u,0,0) + val(u,1,0);
    double vc = val(v,0,0) + val(v,0,1);
    val(ke,0,0) = (uc*uc + vc*vc)/8.;
#else
    double uc = val(u,0,0)*val(u,0,0) + val(u,1,0)*val(u,1,0);
    double vc = val(v,0,0)*val(v,0,0) + val(v,0,1)*val(v,0,1);
    val(ke,0,0) = (uc + vc)/4.;
#endif
    val(psi,0,0) = (val(v,0,0) - val(v,-1,0) + val(u,0,-1) - val(u,0,0))/DX;
  }
}

void advance (void * grid, double t, var * f, var * df)
{
  var u = f[0], v = f[1], h = f[2];
  var du = df[0], dv = df[1], dh = df[2];

  advection_centered (grid, h, u, v, dh);
  momentum (grid, u, v, h, du, dv);
}

#include "init.h"

void update (void * grid, double t, var * f)
{
  var u = f[0], v = f[1], h = f[2];
  boundary_h (grid, h);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);
  boundary_ke_psi (grid);
}

int main (int argc, char ** argv)
{
  #include "parameters.h"

  double t = 0;
  int i = 0;

  void * grid = init_grid (N);
  initial_conditions (grid);
  boundary_b (grid);
  boundary_h (grid, h);
  boundary_u (grid, u, v);
  ke_psi (grid, u, v);
  boundary_ke_psi (grid);

  clock_t start, end;
  start = clock ();
  do {
    double dt = timestep (grid);
    #include "output.h"
#if 1
    advection_centered (grid, h, u, v, hn);
    foreach (grid) { val(h,0,0) += val(hn,0,0)*dt; }
    boundary_h (grid, h);
    momentum (grid, u, v, h, un, vn);
    foreach (grid) {
      val(u,0,0) += val(un,0,0)*dt;
      val(v,0,0) += val(vn,0,0)*dt;
    }
    boundary_u (grid, u, v);
    ke_psi (grid, u, v);
    boundary_ke_psi (grid);
#else /* unstable! */
    var f[3] =  { u, v, h };
    var df[2][3] = {{ un,  vn,  hn },
		    { var(un1), var(vn1), var(hn1) }};
    runge_kutta (2, grid, t, dt, 3, f, df, advance, update);
#endif
    t += dt; i++;
  } while (t < TMAX && i < IMAX);
  end = clock ();
  double cpu = ((double) (end - start))/CLOCKS_PER_SEC;
  fprintf (stderr, "# " GRIDNAME ", %d timesteps, %g CPU, %.3g points.steps/s\n",
	   i, cpu, (N*N*(double)i/cpu));

  free_grid (grid);
}
