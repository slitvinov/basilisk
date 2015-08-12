int coarsen_func (Point point)
{
  return level > 5;
}

int main()
{
  init_grid (64);
  scalar a[];
  fprintf (stderr, "depth: %d mem: %ld\n", depth(), pmtrace.total);
  refine (level < 8, NULL);
  fprintf (stderr, "depth: %d mem: %ld\n", depth(), pmtrace.total);
  coarsen_function (coarsen_func, NULL);
  long tnc = 0;
  foreach()
    tnc++;
  fprintf (stderr, "depth: %d leaves: %ld mem: %ld\n",
	   depth(), tnc, pmtrace.total);
}
