static int ilast = -1;
static double tlast = -1;
static timer progresst;

static void last_events()
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  t = 0., iter = 0;
  double dt = 1e-5;
  while (events (false) && t < HUGE) {
    dtnext (dt);
    dt *= 2.;
    iter = inext, t = tnext;
  }
  ilast = iter, tlast = t;
  if (tlast > 1e10)
    tlast = 0.;
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
  progresst = timer_start();
}

event progress (i += 5)
{
  static FILE * fp = fopen ("progress", "w");
  double peri = ilast ? i/(double)ilast : 1.,
    pert = tlast ? t/tlast : 1., per = min(peri, pert);
  double rem = per ? timer_elapsed (progresst)*(1. - per)/per : 0.;
  fprintf (fp, "%2.0f%% done", per*100.);
  if (rem > 0.) {
    int min = rem/60;
    fprintf (fp, ", %02d:%02d remaining\r", min, (int)rem - 60*min);
  }
  else
    fputc ('\r', fp);
  fflush (fp);
}
