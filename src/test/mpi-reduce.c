#include "utils.h"
#define arr_size 10

int main ()
{
  init_grid (64);

  scalar s[];
  foreach()
    s[] = x + y;
  boundary ({s});

  // statsf() uses reduction operations
  stats stat = statsf (s);
  fprintf (qerr, "%g %g %g\n", stat.min, stat.sum, stat.max);

  // Array reduction
  int cells[arr_size] = {0};
  foreach (reduction(+:cells[:arr_size])) 
    cells[(int)(10*fabs(x))]++;

  for (int i = 0; i < arr_size; i++) 
    fprintf (qerr, "%d ", cells[i]);
  fputc ('\n', qerr);
}
