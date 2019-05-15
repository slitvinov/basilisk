static int _ilast = -1;
static double _tlast = -1;
static timer _progresst;

static void last_events()
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  t = 0., iter = 0;
  // look for events happening at specific times
  while (events (false) && tnext < HUGE)
    t = tnext, iter = inext;
  if (t > 0.)
    _ilast = 0, _tlast = t;
  else {
    // look for condition on maximum time
    double tmin = 0., tmax = HUGE;
    while ((tmax - tmin)/(tmax + tmin) > 0.01) {
      t = (tmin + tmax)/2.;
      if (events (false))
	tmin = t;
      else
	tmax = t;
    }
    if (t < HUGE/2.)
      _ilast = 0, _tlast = t;
    else {
      // look for condition on maximum number of iterations
      _tlast = t = 0.;
      iter = 0;
      while (events (false) && iter < 1<<16)
	iter = inext;
      _ilast = iter < 1<<16 ? iter : 0;
    }
  }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
  _progresst = timer_start();
#if 0
  fprintf (stderr, "_tlast: %g _ilast: %d\n", _tlast, _ilast);
  exit (1);
#endif
}

event progress (i += 5)
{
  static FILE * fp = fopen ("progress", "w");
  double peri = _ilast ? i/(double)_ilast : 1.,
    pert = _tlast ? t/_tlast : 1., per = min(peri, pert);
  if (per > 0.01 && per < 1.) {
    fprintf (fp, "%2.0f%% done", floor(per*100.));
    double rem = per ? timer_elapsed (_progresst)*(1. - per)/per : 0.;
    if (rem > 0.) {
      double min = floor(rem/60.);
      fprintf (fp, ", %02.0f:%02.0f remaining\r", min, rem - 60.*min);
    }
    else
      fputc ('\r', fp);    
  }
  else
    fputc ('\r', fp);
  fflush (fp);
}
