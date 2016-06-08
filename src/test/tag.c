#include "utils.h"
#include "tag.h"

int main()
{
  init_grid (256);
  X0 = Y0 = -0.5;
  scalar t[];
  foreach()
    t[] = sin(4*pi*x)*cos(4.*pi*y) > 0.5;
  tag (t);
  output_ppm (t, file = "t.ppm");
}
