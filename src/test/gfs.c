#include "utils.h"
double t = 0;
#include "input.h"

int main()
{
  init_grid (32);
  scalar a[];
  vector u[];
  foreach()
    u.x[] = u.y[] = a[] = sin(2.*pi*x)*cos(2.*pi*y);

  FILE * fp = fopen("test.gfs", "w");
  output_gfs (fp);
  fclose (fp);

  fp = fopen("test.gfs", "r");
  input_gfs (fp);
  fclose (fp);

  foreach() {
    assert (fabs(a[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
    assert (fabs(u.x[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
    assert (fabs(u.y[] - sin(2.*pi*x)*cos(2.*pi*y)) < 1e-12);
  }    
}
