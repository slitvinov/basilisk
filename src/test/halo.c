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

  foreach_halo()
    fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
  foreach_halo_fine_to_coarse()
    fprintf (stderr, "%g %g %d %d res\n", x, y, level, cell.neighbors);
  output_cells (stdout);

  free_grid ();
}
