// reads the output of gfs.c

double t = 0;
#include "input.h"

int main()
{
  scalar a[];
  vector u[];
  
  FILE * fp = fopen("../gfs/test.gfs", "r");
  input_gfs (fp);
  fclose (fp);

  output_cells (stdout);
  
  foreach() {
    assert (fabs(a[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
    assert (fabs(u.x[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
    assert (fabs(u.y[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
  }

  system ("gfsview-batch2D ../gfs/test.gfs < /dev/null");
}
