/**
# Large-amplitude standing wave

This test case is discussed in [Popinet
(2019)](/Bibliography#popinet2019) and compares the results obtained
with the layered model with the Navier-Stokes/VOF solver for the
large-amplitude, non-linear oscillation of a gravity wave.

~~~gnuplot Free surface profiles
set xlabel 'x'
set ylabel 'z'
plot 'log' u 1:2 w l t 'multilayer', \
     '../large-ns/log' u 1:2 w l t 'VOF', \
     'large.nometric' u 1:2 w l t 'multilayer (no metric)'
~~~

~~~gnuplot Maximum vertical velocity and slope
set xlabel 't'
set ylabel ''
set yr [:0.8]
plot 'out' u 1:2 w l t 'Maximum vertical velocity', \
     'out' u 1:3 w l t 'Maximum slope', \
     '../large-ns/out' u 1:2 every 10 t '', \
     '../large-ns/out' u 1:3 every 10 t '', \
     9.81*x t '9.81 t'
~~~

## See also

* [Same test with VOF](large-ns.c)
*/

#include "grid/multigrid1D.h"
#include "layered/hydro1.h"
#include "layered/nh-box1.h"
#include "layered/remap.h"

int main()
{
  N = 512;
  G = 9.81;
  nl = 20;
  run();
}

event init (i = 0)
{
  foreach() {
    zb[] = - 0.5;
    int l = 0;
    for (scalar h in hl)
      h[] = beta[l++]*(0.07*cos(2.*pi*x) - zb[]);
  }
}

#if 0
event gnuplot (i += 1) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11; set size ratio -1\n");
  fprintf (fp,
	   "set title 'nl = %d, t = %.2f'\n"
	   "p [%g:%g][-0.5:]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 t '' w filledcu lc -1", nl, t,
	   X0, X0 + L0);
  int i = 4;
  for (scalar h in hl)
    fprintf (fp, ", '' u 1:%d w l lw 2 t ''", i++);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    for (scalar h in hl)
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    double z = zb[];
    for (scalar h in hl) {
      fprintf (fp, " %g", z);
      z += h[];
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 0.01\n");
  fflush (fp);
}
#endif

event logfile (i++) {
  double wmax = 0., smax = 0.;
  foreach (reduction (max:wmax) reduction (max:smax)) {
    scalar w, h;
    double Hm = 0., Hp = 0.;
    for (w,h in wl,hl) {
      if (w[] > wmax)
	wmax = w[];
      Hm += h[-1], Hp += h[1];
    }
    if ((Hp - Hm)/(2.*Delta) > smax)
      smax = (Hp - Hm)/(2.*Delta);
  }
  printf ("%g %g %g\n", t, wmax, smax);
}

event profiles (t = 0.1; t += 0.1; t <= 0.5)
{
  foreach() {
    double eta = zb[];
    for (scalar h in hl)
      eta += h[];
    fprintf (stderr, "%g %g\n", x, eta);
  }
  fprintf (stderr, "\n");
}
