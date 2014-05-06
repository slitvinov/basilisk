// Catching floating-point-exceptions

@include <signal.h>

static void caught_abort (int sig)
{
  fprintf (stderr, "Caught signal %d (Aborted)\n", sig);
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', stderr);
  }
}

static void caught_fpe (int sig)
{
  fprintf (stderr, "Caught signal %d (Floating Point Exception)\n", sig);
  abort();
}

static void caught_segfault (int sig)
{
  fprintf (stderr, "Caught signal %d (Segmentation Fault)\n", sig);
  abort();
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL); // FPE
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL); // Segfault
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL); // Abort
}
