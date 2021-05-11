#include "layered/check_eta.h"

event gnuplots (i += 10)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n");
  FILE * fp1 = fopen ("gnuplot", "w");
  foreach()
    fprintf (fp1, "%g %g %g %g\n", x, eta[], deta[], res_eta[]);
  fclose (fp1);
  fprintf (fp,
	   "p 'gnuplot' u 1:3 w l t 'etap - eta',"
	   "  '' u 1:4 w p t 'res'\n");
  fflush (fp);
}
