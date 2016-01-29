#include "saint-venant.h"

void check_flag()
{
  int flag = 1 << user, flag1 = 1 << (user + 1), flag2 = 1 << (user + 2);
  bool failed = false;
  foreach_cell_post (!is_leaf(cell)) {
    if (cell.flags & flag1) {
      fprintf (stderr, "fail %g %g\n", x, y);
      failed = true;
    }
    assert (!(cell.flags & flag));
    //    assert (!(cell.flags & flag1));
    assert (!(cell.flags & flag2));
    assert (!(cell.flags & (1 << (user + 3))));
  }
  if (failed) {
    output_cells (stdout);
    assert (false);
  }
}

@if _MPI
#include "grid/balance.h"
@endif

int LEVEL = 7;

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

event outputfile (t <= 2.5; t += 2.5/8) {
#if !_MPI
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({eta}, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf++);
  output_field ({l});

  /* check symmetry */
  foreach() {
    double h0 = h[];
    point = locate (-x, -y);
    //    printf ("%g %g %g %g %g\n", x, y, h0, h[], h0 - h[]);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (-x, y);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (x, -y);
    assert (fabs(h0 - h[]) < 1e-12);
  }
#endif
}

#if _MPI
event image(i++)
{
  scalar pid[];
  foreach()
    pid[] = pid();
  static FILE * fp = fopen ("pid", "w");
  output_ppm (pid, fp, min = 0, max = npe() - 1);
}
#endif

event adapt (i++) {
  check_flag();
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  check_flag();
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
@if _MPI
  if (s.nf || s.nc) {
    do ; while(balance (0));
    boundary (all);
  }
@endif
}

#if 0
event neighborhood2 (i++) {
  scalar index[];
  indexing (index, 0.);
  char name[200];
  sprintf (name, "neigh-%d", pid());
  FILE * fp = fopen (name, "w");
  neighborhood1 (index, pid(), fp);
  fclose (fp);
#if 1
  sprintf (name,
	   "awk '{print $1,$2}' mpi-restriction-rcv-%d | sort > res-%d;"
	   "awk '{print $1,$2}' neigh-%d | sort | uniq > nei-%d;"
	   "diff nei-%d res-%d > diff-%d; cat diff-%d >> /dev/stderr;",
	   pid(), pid(), pid(), pid(), pid(), pid(), pid(), pid());
  system (name);
#endif
}
#endif
