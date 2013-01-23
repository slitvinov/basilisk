#include <math.h>
#include <stdio.h>
#include <time.h>

struct _Data {
  double a, b;
};

#include GRID

int main ()
{
  for (int l = 6; l <= 11; l++) {
    void * grid = init_grid (1 << l);
    var a = var(a), b = var(b);
    int nloops, i;
    clock_t start, end;

    nloops = i = (1 << 25) >> 2*l;
    start = clock();
    while (i--)
      foreach(grid)
	val(b,0,0) = (val(a,0,1) + val(a,1,0) + val(a,0,-1) + val(a,-1,0) - 4.*val(a,0,0));
    end = clock();
    fprintf (stderr, "%d %g\n", l, 1e9*(end - start)/(double)CLOCKS_PER_SEC/(nloops*(1 << 2*l)));

    free_grid(grid);
  }
}
