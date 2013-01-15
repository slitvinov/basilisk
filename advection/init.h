#include <stdio.h>

#define R0 0.1

void initial_conditions (Data * m, int n)
{
  foreach (m) {
    double x = XC(I), y = YC(J);
    hn(0,0) = exp(-(x*x + y*y)/(R0*R0));
    un(0,0) = 1.;
    vn(0,0) = 0.333;
  }
}

void output_field (Data * m, int n, FILE * fp)
{
  fprintf (fp, "# 1:x 2:y 3:F\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g\n", XC(i), YC(j), M(i,j).h);
    fprintf (fp, "\n");
  }
}
