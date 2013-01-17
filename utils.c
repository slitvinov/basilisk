#include <math.h>

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
// end iterations
int IMAX = 1 << 30;
// CFL number
double CFL = 0.5;

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

void symmetry (Data * m, int n, var v)
{
  foreach_boundary (m, n, right)  { stencil(v,+1,0) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, left)   { stencil(v,-1,0) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, top)    { stencil(v,0,+1) = val(v); } end_foreach_boundary();
  foreach_boundary (m, n, bottom) { stencil(v,0,-1) = val(v); } end_foreach_boundary();  
}

void uv_symmetry (Data * m, int n, var u, var v)
{
  foreach_boundary (m, n, right) {
    stencil(u,+1,0) = 0.;
    stencil(v,+1,0) = val(v);
  } end_foreach_boundary();
  foreach_boundary (m, n, left) {
    stencil(u,-1,0) = 0.;
    stencil(v,-1,0) = val(v);
  } end_foreach_boundary();
  foreach_boundary (m, n, top) {
    stencil(v,0,+1) = 0.;
    stencil(u,0,+1) = val(u);
  } end_foreach_boundary();
  foreach_boundary (m, n, bottom) {
    stencil(v,0,-1) = 0.;
    stencil(u,0,-1) = val(u);
  } end_foreach_boundary();
}

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

