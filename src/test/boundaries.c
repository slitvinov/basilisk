static int refine_func (Point point, void * data)
{
  return (x*x + y*y < 0.25*0.25 ||
	  sq(x - 0.5) + sq(y - 0.5) < 0.25*0.25);
}

int main ()
{
  X0 = -0.5;
  init_grid(8);
  refine_function (refine_func, NULL, NULL);
  output_cells (stdout);
  foreach_boundary_face (bottom)
    fprintf (stderr, "bottom %g %g %d\n", x - Delta/2., y, level);
  foreach_boundary_face (right)
    fprintf (stderr, "right %g %g %d\n", x, y - Delta/2., level);
  foreach_boundary_cell_post (right, !is_leaf (cell))
    fprintf (stderr, "post %g %g %d\n", x, y, level);
  foreach_boundary_halo (right)
    fprintf (stderr, "halo %g %g %d\n", x, y, level);
  free_grid();
}
