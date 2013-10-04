static int refine_func (Point point, void * data)
{
  return x*x + y*y < 0.25*0.25;
}

int main()
{
  X0 = Y0 = -0.5;
  init_grid(8);
  refine_function (refine_func, NULL, NULL);
  output_cells (stdout);
  foreach_vertex()
    fprintf (stderr, "%g %g %d\n", x, y, level);
  free_grid();
}
