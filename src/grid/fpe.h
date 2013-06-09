// Catching floating-point-exceptions

#include <signal.h>

static Point last_point;

static void caught_fpe (int sig)
{
  fprintf (stderr, "Caught signal %d (Floating Point Exception)\n", sig);
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', stderr);
  }
  abort();
}

void catch_fpe (void)
{
  struct sigaction act = {};
  act.sa_handler = caught_fpe;
  last_point.level = -1;
  sigaction (8, &act, NULL);
}
