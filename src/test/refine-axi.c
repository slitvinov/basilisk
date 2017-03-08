// tests that axisymmetric metric is correctly refined

#include "axi.h"
#include "run.h"

double R0 = 0.1;
int LEVEL = 6;

int main()
{
  L0 = 2;
  N = 16;
  run();
}

event init (t = 0) {
  refine (level < LEVEL && sq(x - 2.) + sq(y) < sq(1.5*R0));

  output_cells (stdout);
  
  foreach()
    assert (fabs (cm[] - y) <= 1e-20);
  
  foreach_face(y) {    
    // fprintf (stderr, "%g %g %g %g\n", x, y, fm.y[], fm.y[] - y);
    assert (fabs (fm.y[] - y) <= 1e-20);
  }
  
  foreach_face(x) {
    // fprintf (stderr, "%g %g %g %g\n", x, y, fm.x[], fm.x[] - y);
    assert (fabs (fm.x[] - y) <= 1e-20);
  }
}
