#include "grid/cartesian1D.h"
#include "conservation.h"

scalar h[], q[];
scalar * conserved = {h, q};

q[left]  = - q[];
q[right] = - q[];

void parameters()
{
  theta = 1.;
  X0 = Y0 = -0.5;
  N = 500;
}

double G = 1.;

void flux (const double * s, double * f, double e[2])
{
  double h = s[0], q = s[1], u = q/h;
  f[0] = q;
  f[1] = q*u + G*h*h/2.;
  // min/max eigenvalues
  double c = sqrt(G*h);
  e[0] = u - c; // min
  e[1] = u + c; // max
}

void init()
{
  foreach()
    h[] = 0.1 + exp(-200*x*x);
}

event logfile (t += 0.1; t <= 0.7) {
  foreach()
    fprintf (stderr, "%g %g %g\n", x, h[], q[]);
  fprintf (stderr, "\n");
}

int main() { run(); }
