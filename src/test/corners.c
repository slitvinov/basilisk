int main ()
{
  init_grid(2);
  scalar s[];
  foreach()
    s[] = 1.;
  boundary ({s});
  restriction ({s});
  boundary_restriction ({s});

  for (int l = 0; l <= 1; l++) {
    fprintf (stderr, "level %d\n", l);
    foreach_level(l)
      output_stencil (s, stderr);
  }

  free_grid();
}
