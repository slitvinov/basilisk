int py_scalar_init (scalar s, PyObject * f) {
  if (!PyCallable_Check(f)) {
    fprintf (stderr, "parameter must be callable\n");
    return -1;
  }
  PyObject * code = PyObject_GetAttrString (f, "func_code");
  PyObject * nargs = PyObject_GetAttrString (code, "co_argcount");
  Py_DECREF (code);
  int n = PyInt_AsLong (nargs);
  Py_DECREF (nargs);
  char * format = (n == 0 ? "()" :
		   n == 1 ? "(d)" : 
		   n == 2 ? "(dd)" : 
		   n == 3 ? "(ddd)" :
		   "(error)");
  foreach() {
    PyObject * arglist = Py_BuildValue (format, x, y); // z
    PyObject * result = PyEval_CallObject (f, arglist);
    Py_DECREF (arglist);
    if (result == NULL)
      return -1;
    s[] = PyFloat_AsDouble (result);
    Py_DECREF (result);
  }
  boundary ({s});
  return 0;
}

typedef struct {
  PyObject * i, * t, * action;
} PyEvent;

static double get_double (PyObject * item) {
  if (PyFloat_Check (item))
    return PyFloat_AsDouble (item);
  else if (PyInt_Check (item))
    return PyInt_AsLong (item);
  else {
    fprintf (stderr, "expecting a float\n");
    exit (1);
  }
  return 0;
}

static int py_start (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  if (p->i != Py_None) {
    if (PyInt_Check (p->i))
      *i = PyInt_AsLong (p->i);
    else {
      PyObject * item = PySequence_GetItem (p->i, (ev->a = 0));
      if (!item || !PyInt_Check(item)) {
	fprintf (stderr, "expecting an integer\n");
	exit (1);
      }
      *i = PyInt_AsLong (item);
      Py_DECREF (item);
    }
  }
  else {
    if (PyFloat_Check (p->t))
      *t = PyFloat_AsDouble (p->t);
    else if (PyInt_Check (p->t))
      *t = PyInt_AsLong (p->t);
    else {
      PyObject * item = PySequence_GetItem (p->t, (ev->a = 0));
      if (!item) {
	fprintf (stderr, "expecting a float\n");
	exit (1);
      }
      *t = get_double (item);
      Py_DECREF (item);
    }    
  }
  return 0;
}

static int py_end (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  if (p->i != Py_None)
    return ev->a < PySequence_Size (p->i);
  else
    return ev->a < PySequence_Size (p->t);
}

static int py_inc (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  ev->a++;
  if (p->i != Py_None) {
    if (ev->a < PySequence_Size (p->i)) {
      PyObject * item = PySequence_GetItem (p->i, ev->a);
      *i = PyInt_AsLong (item);
      Py_DECREF (item);
    }
    else
      (*i)++;
  }
  else {
    if (ev->a < PySequence_Size (p->t)) {
      PyObject * item = PySequence_GetItem (p->t, ev->a);
      *t = get_double (item);
      Py_DECREF (item);
    }
    else
      (*t)++;
  }
  return 0;
}

static int py_action (const int i, const double t, Event * ev) {
  PyEvent * p = ev->data;
  PyObject * arglist = Py_BuildValue ("(id)", i, t);
  PyObject * result = PyEval_CallObject (p->action, arglist);
  Py_DECREF (arglist);
  if (result == NULL)
    return 0;
  Py_DECREF (result);
  return 0;
}

int py_register_event (PyObject * action, PyObject * i, PyObject * t) {
  if (!PyCallable_Check(action)) {
    fprintf (stderr, "parameter must be callable\n");
    return -1;
  }

  Event ev;
  ev.last = false;
  ev.action = py_action;
  ev.arrayi = NULL;
  ev.arrayt = NULL;
  PyObject * code = PyObject_GetAttrString (action, "__code__");
  PyObject * file = PyObject_GetAttrString (code, "co_filename");
  ev.file = strdup (PyString_AsString (file));
  Py_DECREF (file);
  PyObject * lineno = PyObject_GetAttrString (code, "co_firstlineno");
  ev.line = PyInt_AsLong (lineno);
  Py_DECREF (lineno);
  PyObject * name = PyObject_GetAttrString (action, "func_name");
  ev.name = strdup (PyString_AsString (name));
  Py_DECREF (name);
  ev.a = 0;

  if (i != Py_None && !PySequence_Check(i) && !PyInt_Check(i)) {
    fprintf (stderr, "parameter i must be a sequence or an int\n");
    exit (1);
  }
  if (t != Py_None && !PySequence_Check(t) && !PyFloat_Check(t)) {
    fprintf (stderr, "parameter t must be a sequence or a float\n");
    exit (1);
  }

  if (PyInt_Check(i) || PyFloat_Check(t)) {
    ev.expr[0] = py_start;
    ev.nexpr = 1;
  }
  else { // sequence
    ev.expr[0] = py_end; ev.expr[1] = py_start; ev.expr[2] = py_inc;
    ev.nexpr = 3;
  }

  // fixme: memory should be freed before exit
  PyEvent * p = ev.data = malloc (sizeof (PyEvent));
  p->i = i; Py_INCREF (i);
  p->t = t; Py_INCREF (t);
  p->action = action; Py_INCREF (action);

  event_register (ev);
  return 0;
}
