int main ()
{
  init_grid(2);
  scalar s = new scalar;
  foreach()
    s[] = 1.;
  restriction (s);

  for (int l = 0; l <= 1; l++)
    boundary_level (scalars (s), l);

  for (int l = 0; l <= 1; l++) {
    fprintf (stderr, "level %d\n", l);
    foreach_level(l)
      output_stencil (s, stderr);
  }

  free_grid();
}