#define INIT ev->expr[0]
#define COND ev->expr[1]
#define INC  ev->expr[2]
#define END_EVENT 1234567890

static void event_error (Event * ev, const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else 
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
	int i = -123456; double t = i;
	(* ev->expr[j]) (&i, &t, ev);
	if (i == -123456 && t == -123456) {
	  /* nothing done to i and t: this must be the condition */
	  if (cond)
	    event_error (ev, "events can only use a single condition");
	  cond = ev->expr[j];
	}
	else {
	  /* this is either an initialisation or an increment */
	  int i1 = i; double t1 = t;
	  (* ev->expr[j]) (&i1, &t1, ev);
	  if (i1 == i && t1 == t) {
	    /* applying twice does not change anything: this is an
	       initialisation */
	    if (init)
	      event_error (ev, "events can only use a single initialisation");
	    init = ev->expr[j];
	  }
	  else {
	    /* this is the increment */
	    if (inc)
	      event_error (ev, "events can only use a single increment");
	    inc = ev->expr[j];
	  }
	}
      }
      INIT = init;
      COND = cond;
      INC  = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (INIT) {
      (* INIT) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == END_EVENT) {
	ev->i = END_EVENT; ev->t = -1;
      }
    }
    else if (INC) {
      (* INC) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
	ev->i = 0;
      if (ev->t != -1)
	ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void init_events (void) {
  for (Event * ev = Events; !ev->last; ev++) 
    if (!strcmp (ev->name, "defaults")) {
#if DEBUG_EVENTS
      char * root = strstr (ev->file, BASILISK);
      fprintf (stderr, "  %-25s %s%s:%d\n", ev->name, 
	       root ? "src" : "",
	       root ? &ev->file[strlen(BASILISK)] : ev->file, 
	       ev->line);
#endif
      ev->action (0, 0., ev);
      event_finished (ev);
    }
    else
      init_event (ev);
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, found = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (found < 0 && !strcmp (event.name, ev->name))
      found = n;
    n++;
  }
  Events = realloc (Events, (n + 2)*sizeof (Event));
  Events[n + 1].last = true;
  if (found >= 0)
    for (; n > found; n--)
      Events[n] = Events[n-1];
  Events[n] = event;
  init_event (&Events[n]);
}

static int event_cond (Event * ev, int i, double t)
{
  if (!COND)
    return true;
  return (* COND) (&i, &t, ev);
}

static int event_do (Event * ev, int i, double t, bool action)
{
  if ((i > ev->i && t > ev->t) || !event_cond (ev, i, t))
    return event_finished (ev);
  if (i == ev->i || fabs (t - ev->t) <= 1e-9) {
#if DEBUG_EVENTS
    if (action) {
      char * root = strstr (ev->file, BASILISK);
      fprintf (stderr, "  %-25s %s%s:%d\n", ev->name, 
	       root ? "src" : "",
	       root ? &ev->file[strlen(BASILISK)] : ev->file, 
	       ev->line);
    }
#endif
    if (action && (* ev->action) (i, t, ev)) {
      event_finished (ev);
      return event_stop;
    }
    if (ev->arrayi) { /* i = {...} */
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
	return event_finished (ev);
    }
    if (ev->arrayt) { /* t = {...} */
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
	return event_finished (ev);
    }
    else if (INC) {
      int i0 = ev->i;
      (* INC) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
	ev->i += i + 1;
      if (!event_cond (ev, i + 1, ev->t))
	return event_finished (ev);
    }
    return event_alive;
  }
  return event_alive;
}

static void end_event_do (int i, double t, bool action)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      ev->action (i, t, ev);
}

int events (int i, double t, bool action)
{
#if DEBUG_EVENTS
  if (1/*action*/)
    fprintf (stderr, "\nevents (i = %d, t = %g)\n", i, t);
#endif
  int inext = 0, cond = 0, cond1 = 0;
  tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT && 
	(COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, i, t, action);
    if (status == event_stop) {
      end_event_do (i, t, action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
	(COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > i)
      inext = 1;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext))
    return 1;
  end_event_do (i, t, action);
  return 0;
}

double dtnext (double t, double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX); // check that dt is not too small
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
	dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
	dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
