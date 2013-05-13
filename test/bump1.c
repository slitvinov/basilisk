#include "grid/cartesian1D.h"
#include "conservation.h"

scalar h = new scalar, q = new scalar;

// fixme: all this should be automatic
scalar h1 = new scalar, q1 = new scalar;
vector gh = new vector, gq = new vector;
scalar dh = new scalar, dq = new scalar;

q[left]  = -q[];
q[right] = -q[];

void parameters()
{
  theta = 1.;
  X0 = Y0 = -0.5;
  N = 500;

  conserved = list_concat (none, {h, q});

  // fixme: all this should be automatic
  conserved1 = list_concat (none, {h1, q1});
  slopes = malloc (3*sizeof (vector));
  slopes[0] = gh; slopes[1] = gq; slopes[2].x = slopes[2].y = -1;
  tendencies = list_concat (none, {dh, dq});
}

double riemann (state * c, double delta, double * f, double dtmax)
{
  double G = 1.;
  double hm = c[0].r, hp = c[0].l;
  double qm = c[1].r, qp = c[1].l;
  double um = qm/hm, up = qp/hp;
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double a = max(ap, -am);
  if (a > 0.) {
    f[0] = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    f[1] = (ap*(qm*um + G*sq(hm)/2.) - am*(qp*up + G*sq(hp)/2.) + 
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*delta/a;
    if (dt < dtmax)
      dtmax = dt;
  }
  else
    f[0] = f[1] = 0.;
  return dtmax;
}

void init()
{
  foreach()
    h[] = 0.1 + exp(-200*x*x);
}

int event (t += 0.1; t <= 0.7) {
  foreach()
    fprintf (stderr, "%g %g %g\n", x, h[], q[]);
  fprintf (stderr, "\n");
}

int main() { run(); }
