#include "grid/cartesian1D.h"
#include "saint-venant.h"

double e1 = 0., e2 = 0., emax = 0.;
int ne = 0;

double a = 3000.;
double h0 = 10.;
double tau = 1e-3;
double Bf = 5.;

int main()
{
  origin (-5000, 0);
  size (10000.);
  G = 9.81;
  for (N = 32; N <= 512; N *= 2) {
    e1 = e2 = emax = 0.;
    ne = 0;
    run();
    fprintf (stderr, "%d %g %g %g\n", N, e1/ne/h0, sqrt(e2/ne)/h0, emax/h0);
  }
}

double Psi (double x, double t)
{
  // Analytical solution, see Sampson, Easton, Singh, 2006
  double p = sqrt (8.*G*h0)/a;
  double s = sqrt (p*p - tau*tau)/2.;
  return a*a*Bf*Bf*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + 
	    (tau*tau/4. - s*s)*cos (2.*s*t)) - Bf*Bf*exp(-tau*t)/(4.*G) -
    exp (-tau*t/2.)/G*(Bf*s*cos (s*t) + tau*Bf/2.*sin (s*t))*x;
}

event init (i = 0)
{
  foreach() {
    zb[] = h0*(x/a)*(x/a);
    h[] = max(h0 + Psi(x,0) - zb[], 0.);
  }
}

event friction (i++) {
  // linear friction (implicit scheme)
  foreach()
    u.x[] /= 1. + tau*dt;
  boundary ({u.x});
}

scalar e[];

event error (i++) {
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

event field (t = 1500) {
  if (N == 64) {
    foreach()
      printf ("p %g %g %g %g %g\n", x, h[], u.x[], zb[], e[]);
    printf ("p\n");
  }
}

event umean (t += 50; t <= 6000) {
  if (N == 128) {
    double sq = 0., sh = 0.;
    foreach() {
      sq += Delta*h[]*u.x[];
      sh += Delta*h[];
    }
    printf ("s %g %g %f\n", t, sq/sh, sh);
  }
}
