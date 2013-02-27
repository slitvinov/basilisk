#include "events.h"

event (t = 4; t <= 12.23; t += 0.1)
{
  fprintf (stderr, "t: %g i: %d\n", t, i);
}

event (i = 5) fprintf (stderr, "i: %d\n", i);

event (i = 2; i *= 2) fprintf (stderr, "i^2: %d\n", i);

event (i++) fprintf (stderr, "i++: %d\n", i);

event (t += 1; t <= 20) fprintf (stderr, "t++: %g\n", t);

int main (int argc, char * argv[])
{
  events_init ();
  double t = 0.;
  int i = 0;
  while (events (i, t)) {
    t = tnext; i++;
  }
}
