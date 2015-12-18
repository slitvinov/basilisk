/* See ../figures/mpi-refine1.svg */

int main (int argc, char * argv[])
{
  init_grid (16);
  refine (level < 5 && x < 0.315 && y < 0.438 && x > 0.25 && y > 0.375, NULL);
  unrefine (x < 0.25 && y > 0.5, NULL);
  refine (level < 5 && x < 0.315 && y < 0.5 && x > 0.25 && y > 0.438, NULL);

  scalar s[];
  foreach()
    s[] = 0.;
  boundary ({s});
}
