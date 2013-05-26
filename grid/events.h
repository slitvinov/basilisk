#define INIT ev->expr[0]
#define COND ev->expr[1]
#define INC  ev->expr[2]

static int event_cond (Event * ev, int i, double t)
{
  if (!COND)
    return true;
  return (* COND) (&i, &t);
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

static int event_do (Event * ev, int i, double t)
{
  if ((i > ev->i && t > ev->t) || !event_cond (ev, i, t))
    return event_finished (ev);
  if (i == ev->i || fabs (t - ev->t) <= 1e-9) {
    if ((* ev->action) (i, t)) {
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
      (* INC) (&ev->i, &ev->t);
      if (!event_cond (ev, i + 1, ev->t))
	return event_finished (ev);
    }
    return event_alive;
  }
  return event_alive;
}

static void event_error (Event * ev, const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

void init_events (void)
{
  for (Event * ev = Events; !ev->last; ev++) 
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
	  (* ev->expr[j]) (&i, &t);
	  if (i == -123456 && t == -123456) {
	    /* nothing done to i and t: this must be the condition */
	    if (cond)
	      event_error (ev, "events can only use a single condition");
	    cond = ev->expr[j];
	  }
	  else {
	    /* this is either an initialisation or an increment */
	    int i1 = i; double t1 = t;
	    (* ev->expr[j]) (&i1, &t1);
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
      if (INIT)
	(* INIT) (&ev->i, &ev->t);
      else if (INC) {
	(* INC) (&ev->i, &ev->t);
	if (ev->i != -1)
	  ev->i = 0;
	if (ev->t != -1)
	  ev->t = 0;
      }
  }
}

int events (int i, double t)
{
  int inext = 0, cond = 0, cond1 = 0;
  tnext = INFINITY;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt)
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, i, t);
    if (status == event_stop)
      return 0;
    if (status == event_alive &&
	(COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > i)
      inext = 1;
  }
  return (!cond || cond1) && (tnext != INFINITY || inext);
}

double dtnext (double t, double dt)
{
  if (tnext != INFINITY) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX); // check that dt is not too small
    dt = (tnext - t)/(n + 1);
    if (n > 0)
      tnext = t + dt;
  }
  else
    tnext = t + dt;
  return dt;
}
