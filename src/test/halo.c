/* definition of halo cells after coarsening */

scalar h[];

int main (int argc, char ** argv)
{
  init_grid (32);

  origin (-0.5, -0.5);
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
  
  /* initial coarsening */
  
  adapt_wavelet ({h}, (double []){1e-2}, 5);

  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l)
      fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
  for (int l = depth() - 1; l >= 0; l--)
    foreach_halo (restriction, l)
      fprintf (stderr, "%g %g %d %d res\n", x, y, level, cell.neighbors);
  output_cells (stdout);

  free_grid ();
}
