// Conway's game of life: http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life

#include "grid/cartesian.h"
#include "run.h"

scalar a[], b[], age[];

void parameters()
{
  X0 = Y0 = -0.5;
  N = 256;
}

event init (i = 0)
{
  foreach() {
    a[] = (x*x + y*y < sq(0.2))*(noise() > 0.);
    age[] = a[];
  }
  boundary({a});
}

event life (i++)
{
  foreach() {
    int neighbors = 0;
    for (int i = -1; i <= 1; i++)
      for  (int j = -1; j <= 1; j++)
	neighbors += a[i,j];
    neighbors -= a[];
    b[] = a[] ? (neighbors == 2 || neighbors == 3) : (neighbors == 3);
    age[] = b[]*(age[] + 1);
  }
  boundary({b});

  swap (scalar, a, b);
}

event print (i++; i < 1000)
{
  // use 'animate age.ppm' to play this,
  // and use the spacebar to play frame by frame 
  static FILE * fp1 = fopen ("age.ppm", "w");
  scalar m[];
  foreach()
    m[] = age[] ? 1 : -1;
  output_ppm (age, fp1, 256, mask = m, linear = false );
}

int main() { run(); }
