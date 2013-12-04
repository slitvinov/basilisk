/**
# [Conway's game of life](http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life)

We use a cartesian grid and the generic time loop. */

#include "grid/cartesian.h"
#include "run.h"

/**
We need two fields to store the old and new states as well as a field
to store the age of each cell. */

scalar a[], b[], age[];

/**
The lower-left corner is at (-0.5,-0.5) (the default box size is one)
i.e. the domain spans (-0.5,-0.5) (0.5,0.5) and is discretised using
$256^2$ cells. */

void parameters()
{
  X0 = Y0 = -0.5;
  N = 256;
}

/**
We initialise zeros and ones randomly (the `noise()` function returns
random numbers between -1 and 1) in a circle centered on the origin of
radius 2. */

event init (i = 0)
{
  foreach() {
    a[] = (x*x + y*y < sq(0.2))*(noise() > 0.);
    age[] = a[];
  }
  boundary({a});
}

/**
## Animation

We generate images of the age field every 5 timesteps for the first
1000 timesteps of the evolution. */

event movie (i += 5; i < 1000)
{

/**
The filename is `age-000.ppm`, `age-001.ppm` etc... */

  char name[80];
  static int nf = 0;
  sprintf (name, "age-%03d.ppm", nf++);
  FILE * fp = fopen (name, "w");

/**
We mask out dead cells (i.e. cells for which `age` is zero). */

  scalar m[];
  foreach()
    m[] = age[] ? 1 : -1;
  output_ppm (age, fp, 512, mask = m, linear = false );
  fclose (fp);
}

/**
Finally, after completion of the simulation, we create a compressed
GIF animation. */

event gif (t = end) {

/**
We first convert each PPM file into GIF (this is the most
time-consuming operation!)... */

  system ("for f in *.ppm; do"
	  "  convert $f `basename $f .ppm`.gif && rm -f $f;"
	  "done;"

/**
... then use `gifsicle` to create the final compressed animated GIF. */

	  "gifsicle --colors 256 --optimize --delay 1"
	  " --loopcount=0 age-*.gif > age.gif && rm -f age-*.gif");
}

/**
## Game of life algorithm */

event life (i++)
{
  foreach() {

/**
We count the number of live neighbors in a 3x3 neighbourhood. */

    int neighbors = - a[];
    for (int i = -1; i <= 1; i++)
      for  (int j = -1; j <= 1; j++)
	neighbors += a[i,j];

/**
If a cell is alive and surrounded by 2 or 3 neighbors it carries on
living, otherwise it dies. If a cell is dead and surrounded by exactly
3 neighbors it becomes alive. */

    b[] = a[] ? (neighbors == 2 || neighbors == 3) : (neighbors == 3);

/**
The age of live cells is incremented. */

    age[] = b[]*(age[] + 1);
  }
  boundary({b});

/**
Here we swap the old state (`a`) with the new state (`b`). */

  swap (scalar, a, b);
}

/**
The generic `run()` function implements the main time loop. */

int main() { run(); }

/**
![Evolution of the age of cells](life/age.gif) */
