event event1 (t = 4; t <= 12.23; t += 0.1)
{
  fprintf (stderr, "t: %g i: %d\n", t, i);
}

event event2 (i = 5) fprintf (stderr, "i: %d\n", i);

event event3 (i = 2; i *= 2) fprintf (stderr, "i^2: %d\n", i);

event event4 (i++) fprintf (stderr, "i++: %d\n", i);

event event5 (t += 1; t <= 20) fprintf (stderr, "t++: %g\n", t);

int main (int argc, char * argv[])
{
  init_events();
  double t = 0.;
  int i = 0;
  while (events (i, t)) {
    t = tnext; i++;
  }
}
