/* See ../figures/mpi-refine1.svg */

int main (int argc, char * argv[])
{
  init_grid (16);
  refine (level < 5 && y < 0.315 && (1. - x) < 0.438 && y > 0.25 &&
	  (1. - x) > 0.375, NULL);
  unrefine (y < 0.25 && (1. - x) > 0.5, NULL);
  refine (level < 5 && y < 0.315 && (1. - x) < 0.5 && y > 0.25 &&
	  (1. - x) > 0.438, NULL);

  scalar s[];
  foreach()
    s[] = 0.;
  boundary ({s});
}
