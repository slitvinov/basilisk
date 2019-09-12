/**
# Dispersion relations for various models

This test case measures the dependence of the oscillation frequency of
a gravity wave of wavenumber $k$ on the water depth $h$. The
corresponding phase velocity $c=\omega/k$ is then compared to the
exact linear dispersion relation
$$
c_e = \sqrt{gk\tanh(kh)}
$$

We use a 1D grid and either the Green-Naghdi solver or the multi-layer
solver. */

#include "grid/multigrid1D.h"
#if GN
  #include "green-naghdi.h"
#else
  #define NL 1
  #include "layered/hydro.h"
  #include "layered/nh-box.h"
  // necessary for stability at high h0 of the box scheme
  #include "layered/remap.h"
#endif

double emax = 1.;
double h0 = 0.1;

/**
We initialise a wave of small relative amplitude ($10^{-3}$) to be in
the linear regime. For the multilayer model, we also run a case with
three "optimised" non-uniform layer thicknesses. */

event init (i = 0)
{
#if NL
  if (nl == 3 && emax == 1.5)
    beta[0] = 0.68, beta[1] = 0.265, beta[2] = 0.055;
  foreach() {
    int l = 0;
    for (scalar h in hl)
      h[] = beta[l++]*h0*(1. + 0.001*cos(x));
  }
#else
  foreach()
    h[] = h0*(1. + 0.001*cos(x));
#endif
}

/**
The code below locates the time of zero "up-crossing" of the amplitude
in the middle of the wave and computes the corresponding averaged
wave period. */

double Tm = 0.;
int nm = 0;

event logfile (i++) {
  Point point = locate (L0/2.);
#if NL
  double H = 0.;
  for (scalar h in hl)
    H += h[];
  double dh = (H - h0)/h0;
#else
  double dh = (h[] - h0)/h0;
#endif
  //  static FILE * fp = fopen ("timeseries", "w");
  //  fprintf (fp, "%g %g %g\n", t, H, dh);
  static double told = 0., hold = 0., tsold = -1;
  if (i == 0) {
    Tm = 0., nm = 0;
    told = 0., hold = 0., tsold = -1;
  }
  else {
    if (i > 1 && hold < 0. && dh > 0.) {
      // this is an (upward) zero-crossing at time ts
      double ts = told - hold*(t - told)/(dh - hold);
      //      fprintf (fp, "ts %g 0\n", ts);
      if (tsold > 0. && nm < 4) {
	// we average the periods
	Tm += ts - tsold;
	nm++;
      }
      tsold = ts;
    }
    told = t;
    hold = dh;
  }
}

/**
After ten (exact) wave periods, we stop and dump the solution. */

event dumps (t = 10.*2.*pi/sqrt(tanh(h0))) // ten wave periods
{
  FILE * fp = fopen ("prof", "w");
  fprintf (fp, "x");
  for (scalar s in all)
    fprintf (fp, " %s", s.name);
  fprintf (fp, "\n");
  foreach() {
    fprintf (fp, "%g", x);
    for (scalar s in all)
      fprintf (fp, " %g", s[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "\n");
  fclose (fp);
}

#if 0
event movie (i += 1) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fputs ("set term x11\n", fp);
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "unset key\n"
	   "p []'-' u 1:3 w l, '' u 1:5 w l, '' u 1:7 w l, '' u 1:9 w l\n", t);
  foreach() {
    fprintf (fp, "%g", x);
    scalar h;
    vector q;
    for (h, q in hl, ql)
      fprintf (fp, " %g %g", h[], q.x[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);

  scalar eta[];
  foreach() {
    eta[] = zb[];
    for (scalar h in hl)
      eta[] += h[];
  }
  fprintf (stderr, "%g %g\n", t, statsf(eta).max); 
}
#endif

/**
We output the wave period and compute the relative error `emax`. */

event end (t = end)
{
  fprintf (stderr, "%g %g %g\n", h0, sqrt(tanh(h0)), 2.*pi/(Tm/nm));
  fflush (stderr);
  emax = 2.*pi/(Tm/nm)/sqrt(tanh(h0));
}

/**
We run for several layers and many water depths, stopping when the
error becomes too large. */

int main()
{
  periodic (right);
  size (2.*pi);
  N = 128;
#if NL
  for (nl = 1; nl <= 5; nl++) {
    char name[80];
    sprintf (name, "log-%d", nl);
    freopen (name, "w", stderr);
    emax = 1.;
    DT = 2.*pi/sqrt(tanh(h0))/100.;
    for (h0 = 0.1; h0 <= 100. && emax > 0.96; h0 *= 1.3)
      run();
  }
  freopen ("log-opt", "w", stderr);
  nl = 3;
  emax = 1.;
  for (h0 = 0.1; h0 <= 100. && emax > 0.96; h0 *= 1.3) {
    emax = 1.5;
    run();
  }
#else
  alpha_d = 1.159;
  for (h0 = 0.1; h0 <= 100. && emax > 0.96 && emax < 1.04; h0 *= 1.2)
    run();
#endif
}

/**
Note that the standard Green-Naghdi model has a [2,2] Padé approximant
dispersion relation (see e.g. Clamond, Dutykh, Mitsotakis, 2016). The
same as Nwogu (1993). Here we used the optimized version of Chazel et
al. (2011).

~~~gnuplot Dispersion relation
g = 1.
k = 1.

omega1_keller(h_0)=sqrt((4*h_0*g*k**2)/(h_0**2*k**2+4))

omega2_keller(h_0,h_1)=sqrt(((4*h_0*h_1**2+4*h_0**2*h_1)*g*k**4+(16*h_1+16*h_0)*g*k**2)/(h_0**2*h_1**2*k**4+(4*h_1**2+16*h_0*h_1+4*h_0**2)*k**2+16))

omega3_keller(h_0,h_1,h_2)=sqrt((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*g*k**6+((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*g*k**4+(64*h_2+64*h_1+64*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*k**6+((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*k**4+(16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*k**2+64))

omega4_keller(h_0,h_1,h_2,h_3)=sqrt(((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*h_3**2+4*h_0**2*h_1**2*h_2**2*h_3)*g*k**8+(((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*h_3**2+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*g*k**6+((64*h_2+64*h_1+64*h_0)*h_3**2+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3+(64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*g*k**4+(256*h_3+256*h_2+256*h_1+256*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*h_3**2*k**8+(((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*h_3**2+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3+4*h_0**2*h_1**2*h_2**2)*k**6+((16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*h_3**2+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*k**4+(64*h_3**2+(256*h_2+256*h_1+256*h_0)*h_3+64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*k**2+256))

omega5_keller(h_0,h_1,h_2,h_3,h_4)=sqrt((((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*h_3**2+4*h_0**2*h_1**2*h_2**2*h_3)*h_4**2+4*h_0**2*h_1**2*h_2**2*h_3**2*h_4)*g*k**10+((((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*h_3**2+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_4**2+(((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3**2+((64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_3+16*h_0**2*h_1**2*h_2**2)*h_4+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3**2+16*h_0**2*h_1**2*h_2**2*h_3)*g*k**8+(((64*h_2+64*h_1+64*h_0)*h_3**2+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3+(64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_4**2+((64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3**2+((256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_3+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_4+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3**2+((64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_3+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*g*k**6+((256*h_3+256*h_2+256*h_1+256*h_0)*h_4**2+(256*h_3**2+(1024*h_2+1024*h_1+1024*h_0)*h_3+256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_4+(256*h_2+256*h_1+256*h_0)*h_3**2+(256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_3+(256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*g*k**4+(1024*h_4+1024*h_3+1024*h_2+1024*h_1+1024*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*h_3**2*h_4**2*k**10+((((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*h_3**2+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3+4*h_0**2*h_1**2*h_2**2)*h_4**2+(((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3**2+16*h_0**2*h_1**2*h_2**2*h_3)*h_4+4*h_0**2*h_1**2*h_2**2*h_3**2)*k**8+(((16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*h_3**2+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_4**2+(((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3**2+((64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_3+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_4+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3**2+((64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_3+16*h_0**2*h_1**2*h_2**2)*k**6+((64*h_3**2+(256*h_2+256*h_1+256*h_0)*h_3+64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_4**2+((256*h_2+256*h_1+256*h_0)*h_3**2+(256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_3+(256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_4+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3**2+((256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_3+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*k**4+(256*h_4**2+(1024*h_3+1024*h_2+1024*h_1+1024*h_0)*h_4+256*h_3**2+(1024*h_2+1024*h_1+1024*h_0)*h_3+256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*k**2+1024))

alpha=1.159
omega_gn(k,h)=sqrt(g*h*k**2*(1.+(alpha-1.)*(k*h)**2/3.)/(1.+alpha*(k*h)**2/3.))

set xlabel 'kh'
set ylabel 'c/c_e'
set key bottom left
set xr [0.1:100]
set yr [0.96:1.04]
set logscale x
set grid
plot  omega1_keller(x)/sqrt(tanh(x)) t '1 layer', \
      'log-1' u ($1):($3/sqrt(tanh($1))) t '' pt 5, \
      omega2_keller(x/2.,x/2.)/sqrt(tanh(x)) t '2 layers', \
      'log-2' u ($1):($3/sqrt(tanh($1))) t '' pt 7, \
      omega3_keller(x/3.,x/3.,x/3.)/sqrt(tanh(x)) t '3 layers', \
      'log-3' u ($1):($3/sqrt(tanh($1))) t '' pt 9, \
      omega4_keller(x/4.,x/4.,x/4.,x/4.)/sqrt(tanh(x)) t '4 layers', \
      'log-4' u ($1):($3/sqrt(tanh($1))) t '' pt 11, \
      omega5_keller(x/5.,x/5.,x/5.,x/5.,x/5.)/sqrt(tanh(x)) t '5 layers', \
      'log-5' u ($1):($3/sqrt(tanh($1))) t '' pt 13, \
      omega3_keller(0.68*x,0.265*x,0.055*x)/sqrt(tanh(x)) \
          t 'Optimised 3 layers' lc 0, \
      'log-opt' u ($1):($3/sqrt(tanh($1))) t '' pt 15 lc 0, \
      omega_gn(1,x)/sqrt(tanh(x)) t 'Green-Naghdi (1 layer)' lc 4, \
      '../dispersion-gn/log' u ($1):($3/sqrt(tanh($1))) t ''
~~~
*/
