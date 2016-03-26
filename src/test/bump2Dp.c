// very close to bump2D.c but with additional checks for MPI parallelism

#include "saint-venant.h"
#include "check_restriction.h"

#define BGHOSTS 2

int LEVEL = 6;

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init (i = 0)
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (ferr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 2.5; t += 2.5/8);

event image(i++)
{
  scalar pid[];
  foreach()
    pid[] = pid();
  static FILE * fp = fopen ("pid", "w");
  output_ppm (pid, fp, min = 0, max = npe() - 1, map = randomap);
  foreach()
    pid[] = level;
  static FILE * fp1 = fopen ("level", "w");
  output_ppm (pid, fp1, min = 0, max = LEVEL);
}

event adapt (i++) {

#if BGHOSTS == 2
  scalar s[];
  face vector u[];
  quadtree_trash ({s,u});
  foreach()
    s[] = 1;
  foreach_face()
    u.x[] = 1;
  boundary ({s,u});
#endif
  
#if 0  
  debug_iteration = i;

  char name[80];
  sprintf (name, "uf-%d-%d", i, pid());
  FILE * fp = fopen(name, "w");
  foreach_cell() {
    fprintf (fp, "%g %g %g\n", x - Delta/2., y, u.x[]);
    if (allocated(1))
      fprintf (fp, "%g %g %g\n", x + Delta/2., y, u.x[1]);
  }
  fclose (fp);
#endif
  
  astats st = adapt_wavelet ({h}, (double[]){1e-2}, LEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", st.nf, st.nc);
  restriction ({zb}); // fixme: why is it necessary with MPI?

#if BGHOSTS == 2
  foreach()
    foreach_neighbor()
    //      if (cell.pid >= 0) // fixme
	assert (s[] == 1);
  check_restriction (s);
#if 1
  foreach_face()
    for (int i = -2; i <= 2; i++)
      //      if (neighbor(0,i).pid >= 0)
	assert (u.x[0,i] == 1);
#endif
#endif
}
