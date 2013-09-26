int main ()
{
  init_grid(4);
  output_cells (stdout);
  foreach_boundary_face (top)
    fprintf (stderr, "top %g %g\n", x - delta/2., y - delta/2.);
  foreach_boundary_face (right)
    fprintf (stderr, "right %g %g\n", x - delta/2., y - delta/2.);
  free_grid();
}
