#define INIT ev->expr[0]
#define COND ev->expr[1]
#define INC  ev->expr[2]

static int event_cond (Event * ev, int i, double t)
{
  if (!COND)
    return true;
  return (* COND) (&i, &t);
}

static void event_finished (Event * ev)
{
  ev->t = ev->i = -1;
}

static int event_do (Event * ev, int i, double t)
{
  if ((i > ev->i && t > ev->t) || !event_cond (ev, i, t)) {
    event_finished (ev);
    return 0;
  }
  if (i == ev->i || t == ev->t) {
    if ((* ev->action) (i, t)) {
      event_finished (ev);
      return 1;
    }
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
    else if (INC) {
      (* INC) (&ev->i, &ev->t);
      if (!event_cond (ev, i + 1, ev->t))
	event_finished (ev);
    }
  }
  return 0;
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
  int inext = 0, cond = 0;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (COND || (INIT && !COND && !INC) || ev->arrayi || ev->arrayt)
      cond = 1;
  tnext = undefined;
  for (Event * ev = Events; !ev->last; ev++) {
    if (event_do (ev, i, t))
      return 0;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > i && (!INC || COND || !cond))
      inext = 1;
  }
  return tnext != undefined || inext;
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
