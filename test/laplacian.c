#include <math.h>
#include <stdio.h>
#include <time.h>

scalar a = new scalar, b = new scalar;

int main ()
{
  for (int l = 4; l <= 11; l++) {
    init_grid (1 << l);
    int nloops, i;
    clock_t start, end;

    nloops = i = (1 << 25) >> 2*l;
    start = clock();
    while (i--)
      foreach()
	b[] = (a[0,1] + a[1,0] + a[0,-1] + a[-1,0] - 4.*a[]);
    end = clock();
    fprintf (stderr, "%d %g\n", l, 
	     1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)));

    nloops = i = (1 << 25) >> 2*l;
    double sum = 0.;
    start = clock();
    while (i--)
      foreach()
	sum += a[];
    end = clock();
    printf ("%d %g %g\n", l, 1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)), sum);
    free_grid();
  }
}
