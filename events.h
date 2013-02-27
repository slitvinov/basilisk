static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t);
}

static void event_finished (Event * ev)
{
  ev->t = ev->i = -1;
}

static void event_do (Event * ev, int i, double t)
{
  if ((i > ev->i && t > ev->t) || !event_cond (ev, i, t)) {
    event_finished (ev);
    return;
  }
  if (i == ev->i || t == ev->t) {
    (* ev->action) (i, t);
    if (ev->arrayi) { /* i = {...} */
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
	event_finished (ev);
    }
    if (ev->arrayt) { /* t = {...} */
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
	event_finished (ev);
    }
    else if (ev->expr[2]) { /* increment */
      (* ev->expr[2]) (&ev->i, &ev->t);
      if (!event_cond (ev, i + 1, ev->t))
	event_finished (ev);
    }
  }
}

void events_init (void)
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
	      fprintf (stderr, "warning: event condition redefined\n");
	    cond = ev->expr[j];
	  }
	  else {
	    /* this is either an initialisation or an increment */
	    int i1 = i; double t1 = t;
	    (* ev->expr[j]) (&i1, &t1);
	    if (i1 == i && t1 == t) {
	      /* applying twice does not change anything: this is an initialisation */
	      if (init)
		fprintf (stderr, "warning: event initialisation redefined\n");
	      init = ev->expr[j];
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
	ev->nexpr = 0;
      }
      ev->i = ev->t = -1;
      if (ev->expr[0])
	(* ev->expr[0]) (&ev->i, &ev->t);
      else if (ev->expr[2]) {
	(* ev->expr[2]) (&ev->i, &ev->t);
	if (ev->i != -1)
	  ev->i = 0;
	if (ev->t != -1)
	  ev->t = 0;
      }
  }
}

int events (int i, double t)
{
  tnext = undefined;
  for (Event * ev = Events; !ev->last; ev++) {
    event_do (ev, i, t);
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
