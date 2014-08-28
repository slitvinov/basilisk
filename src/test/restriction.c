// checks that restriction works for face vectors on a non-trivial quadtree grid

int refine_circle (Point point, void * data)
{
  int depth = *((int *)data);
  x -= 0.1; y -= 0.1;
  return (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));
}

int main(int argc, char * argv[])
{
  int depth = argc > 1 ? atoi(argv[1]) : 4;
  X0 = Y0 = -0.7;
  init_grid (1);
  while (refine_function (refine_circle, &depth, NULL));

  face vector u[];
  foreach_face()
    u.x[] = 1.;
  restriction ((scalar *){u});

  output_cells (stdout);

  for (int l = 0; l <= depth(); l++)
    foreach_level(l) {
      fprintf (stderr, "%g %g %g %d\n", x - Delta/2., y, u.x[], l);
      assert (u.x[] == 1. && u.y[] == 1. && u.x[1,0] == 1. && u.y[0,1] == 1.);
    }
}
