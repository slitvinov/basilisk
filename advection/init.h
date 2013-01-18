#include <stdio.h>

#define R0 0.1

void initial_conditions (Data * m, int n)
{
  foreach (m, n) {
    h(0,0) = exp(-(x*x + y*y)/(R0*R0));
  } end_foreach();
}

void output_field (Data * m, int n, var v, FILE * fp)
{
  Data * _m = m; int _n = n; /* fixme: this is a hack */
  fprintf (fp, "# 1:x 2:y 3:F\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g\n", XC(i), YC(j), stencil(v,0,0));
    fprintf (fp, "\n");
  }
}
