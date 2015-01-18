#include <sys/resource.h>

void mem()
{
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  fprintf (stderr, "%ld kB\n", usage.ru_maxrss);
}

int coarsen_func (Point point)
{
  return true;
}

int main()
{
  mem();
  init_grid (64);
  mem();
  scalar a[];
  mem();
  fprintf (stderr, "depth: %d\n", depth());
  coarsen_function (coarsen_func, NULL);
  long tnc = 0;
  foreach()
    tnc++;
  fprintf (stderr, "depth: %d leaves: %ld\n", depth(), tnc);
  mem();
}
