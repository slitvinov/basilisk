/**
# Green-Naghdi soliton

The [Green-Naghdi](/src/green-naghdi.h) system of equations admits
solitary wave solutions of the form
$$
h(\zeta) = h_1 + (h_2 - h_1){\rm sech}^2
    \left(\frac{\zeta}{2}\sqrt{\frac{3(h_2-h_1)}{h_2h_1^2}}\right)
$$
$$
u(\zeta) = D\left(1-\frac{h_1}{h(\zeta)}\right)
$$
with $\zeta = x - ct$ and the soliton velocity $c^2=gh_2$. */

#include "grid/multigrid1D.h"
#include "green-naghdi.h"

/**
The domain is 700 metres long, the acceleration of gravity is $10$
m/s$^{-2}$. We compute the solution in one dimension for a number of
grid points varying between 128 and 1024. */

int main()
{
  X0 = -200.;
  L0 = 700.;
  G = 10.;
  for (N = 128; N <= 1024; N *= 2)
    run();
}

/**
We follow [Le Métayer et al, 2010](/src/references.bib#lemetayer2010)
(section 6.1.2) for $h_1$ and $h_2$. */

double h1 = 10, h2 = 12.1;

double sech2 (double x) {
  double a = 2./(exp(x) + exp(-x));
  return a*a;
}

double soliton (double x, double t)
{
  double D = sqrt(G*h2), psi = x - D*t;
  return h1 + (h2 - h1)*sech2 (psi/2.*sqrt(3.*(h2 - h1)/(h2*h1*h1)));
}

event init (i = 0)
{
  double D = sqrt(G*h2);
  foreach() {
    h[] = soliton (x, t);
    u[] = D*(1. - h1/h[]);
  }
}

/**
We output the profiles and reference solution at regular intervals. */

event output (t = {0,7.3,14.6,21.9,29.2}) {
  if (N == 256) {
    foreach()
      fprintf (stdout, "%g %g %g %g\n", x, h[], u[], soliton (x, t));
    fprintf (stdout, "\n");
  }
}

/**
We compute the error between the theoretical and numerical solutions
at $t = 29.2$. */

event error (t = end) {
  scalar e[];
  foreach()
    e[] = h[] - soliton (x, t);
  fprintf (stderr, "%d %g\n", N, normf(e).max/(h2 - h1));
}

/**
~~~gnuplot Depth profiles for N = 256.
set grid
set xlabel 'x'
set ylabel 'z'
plot 'out' u 1:4 w l t 'exact', 'out' w l t 'numerical'
~~~

The method has a second-order rate of convergence as expected.

~~~gnuplot Relative error as a function of resolution.
set logscale
set xlabel 'N'
set ylabel 'max|e|/a'
set xtics 128,2,1024
set grid
fit a*x+b 'log' u (log($1)):(log($2)) via a,b
plot [100:1250]'log' u 1:($2) pt 7 t '', \
     exp(b)*x**a t sprintf("%.0f/N^{%4.2f}", exp(b), -a)
~~~ */
