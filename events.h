static int event_cond (Event * event, void * grid, int i, double t)
{
  if (!event->expr[1])
    return true;
  return (* event->expr[1]) (grid, &i, &t);
}

static void event_do (Event * event, void * grid, int i, double t)
{
  if (!event->action)
    return;
  if ((i > event->i && t > event->t) || 
      !event_cond (event, grid, i, t)) {
    event->action = NULL;
    event->t = event->i = -1;
    return;
  }
  if (i == event->i || t == event->t) {
    (* event->action) (grid, i, t);
    if (event->expr[2]) {
      (* event->expr[2]) (grid, &event->i, &event->t);
      if (!event_cond (event, grid, i + 1, event->t)) {
	event->action = NULL;
	event->t = event->i = -1;
      }
    }
  }
}

void events_init (void * grid)
{
  for (Event * ev = Events; !ev->last; ev++) {
    Expr init = NULL, cond = NULL, inc = NULL;
    for (int j = 0; j < ev->nexpr; j++) {
      int i = -123456; double t = i;
      (* ev->expr[j]) (grid, &i, &t);
      if (i == -123456 && t == -123456) {
	/* nothing done to i and t: this must be the condition */
	if (cond)
	  fprintf (stderr, "warning: event condition redefined\n");
	cond = ev->expr[j];
      }
      else {
	/* this is either an initialisation or an increment */
	int i1 = i; double t1 = t;
	(* ev->expr[j]) (grid, &i1, &t1);
	if (i1 == i && t1 == t) {
	  /* applying twice does not change anything: this is an initialisation */
	  if (init)
	    fprintf (stderr, "warning: event initialisation redefined\n");
	  init = ev->expr[j];
	  ev->i = i; ev->t = t;
	}
	else {
	  /* this is the increment */
	  if (inc)
	    fprintf (stderr, "warning: event increment redefined\n");
	  inc = ev->expr[j];
	}
      }
    }
    ev->expr[0] = init;
    ev->expr[1] = cond;
    ev->expr[2] = inc;
  }
}

int events (void * grid, int i, double t)
{
  tnext = undefined;
  for (Event * ev = Events; !ev->last; ev++) {
    event_do (ev, grid, i, t);
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
  }
  return tnext != undefined;
}

double dtnext (double t, double dt)
{
  if (tnext != undefined) {
    int n = (tnext - t)/dt;
    dt = (tnext - t)/(n + 1);
    if (n > 0)
      tnext = t + dt;
  }
  else
    tnext = t + dt;
  return dt;
}
