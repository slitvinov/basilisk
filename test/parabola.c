#include "grid/cartesian1D.h"
#include "saint-venant1.h"

void parameters()
{
  L0 = 10000.;
  G = 9.81;
}

double a = 3000.;
double h0 = 10.;
double tau = 1e-3;
double Bf = 5.;

double Psi (double x, double t)
{
  // Analytical solution, see Sampson, Easton, Singh, 2006
  double p = sqrt (8.*G*h0)/a;
  double s = sqrt (p*p - tau*tau)/2.;
  return a*a*Bf*Bf*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + 
	    (tau*tau/4. - s*s)*cos (2.*s*t)) - Bf*Bf*exp(-tau*t)/(4.*G) -
    exp (-tau*t/2.)/G*(Bf*s*cos (s*t) + tau*Bf/2.*sin (s*t))*x;
}

void init()
{
  foreach() {
    zb[] = h0*(x/a)*(x/a);
    h[] = max(h0 + Psi(x,0) - zb[], 0.);
  }
}

event(i++) {
  // linear friction (implicit scheme)
  foreach()
    q.x[] /= 1. + tau*dt;
  boundary (q.x);
}

scalar e = new scalar;

double e1 = 0., e2 = 0., emax = 0.;
int ne = 0;

event (i++) {
  foreach()
    e[] = h[] - max(h0 + Psi(x,t) - zb[], 0.);
  norm n = normf (e);
  e1 += n.avg;
  e2 += n.rms*n.rms;
  ne++;
  if (n.max > emax)
    emax = n.max;
  printf ("e %g %g %g %g\n", t, n.avg, n.rms, n.max);
}

event (t = 1500) {
  if (N == 64) {
    foreach()
      printf ("p %g %g %g %g %g\n", x, h[], q.x[], zb[], e[]);
    printf ("p\n");
  }
}

event (t += 50; t <= 6000) {
  if (N == 128) {
    double sq = 0., sh = 0.;
    foreach() {
      sq += DX*q.x[];
      sh += DX*h[];
    }
    printf ("s %g %g %f\n", t, sq/sh, sh);
  }
}

int main() {
  for (N = 32; N <= 512; N *= 2) {
    e1 = e2 = emax = 0.;
    ne = 0;
    run();
    fprintf (stderr, "%d %g %g %g\n", N, e1/ne/h0, sqrt(e2/ne)/h0, emax/h0);
  }
}
