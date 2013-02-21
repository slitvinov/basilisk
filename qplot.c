#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

FILE * gnuplot = NULL;

void cleanup (int number)
{
  pclose (gnuplot);
  system ("rm -f .qplot??????");
  exit (1);
}

int main (int argc, char * argv[])
{
  float fn;
  char command[80] = "gnuplot";
  int i;
  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-') {
      strcat (command, " "); strcat (command, argv[i]);
    }
  gnuplot = popen (command, "w");
  if (!gnuplot) {
    perror ("qplot: could not open pipe: ");
    return 1;
  }
  fprintf (gnuplot, "set terminal wxt noraise\n");
  for (i = 1; i < argc; i++)
    if (argv[i][0] != '-')
      fprintf (gnuplot, "load '%s'\n", argv[i]);
  signal (SIGINT, cleanup);
  while (fread (&fn, sizeof(float), 1, stdin) == 1) {
    char fname[] = ".qplotXXXXXX";
    int fd = mkstemp (fname), n = fn, nr;
    if (fd < 0) {
      perror ("qplot: could not create temporary file: ");
      return 1;
    }
    FILE * fp = fdopen (fd, "w");
    fwrite (&fn, sizeof(float), 1, fp);
    float * v = malloc (sizeof (float)*(n + 1));
    if ((nr = fread (v, sizeof (float), n, stdin)) != n) {
      fprintf (stderr, "qplot: expecting %d y-values only got %d\n", n, nr);
      remove (fname);
      return 1;
    }
    fwrite (v, sizeof (float), n, fp);
    int i;
    for (i = 0; i < n; i++) {
      if ((nr = fread (v, sizeof (float), n + 1, stdin)) != n + 1) {
	fprintf (stderr, "qplot: expecting %d x-z-values only got %d\n", n + 1, nr);
	remove (fname);
	return 1;
      }
      fwrite (v, sizeof (float), n + 1, fp);
    }
    free (v);
    fclose (fp);
    fprintf (gnuplot, "splot '%s' binary u 2:1:3 t ''\n", fname);
    fprintf (gnuplot, "! rm -f %s\n", fname);
    fflush (gnuplot);
  }
  pclose (gnuplot);
  system ("rm -f .qplot??????");
  return 0;
}
