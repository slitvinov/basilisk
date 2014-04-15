#include <Python.h>

int py_scalar_init (scalar s, PyObject * f) {
  if (!PyCallable_Check(f)) {
    fprintf (stderr, "parameter must be callable");
    return -1;
  }
  foreach() {
    PyObject * arglist = Py_BuildValue ("(dd)", x, y);
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
  PyObject * i;
  PyObject * action;
} PyEvent;

static int py_start (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  if (PyInt_Check (p->i))
    *i = PyInt_AsLong (p->i);
  else {
    PyObject * item = PyList_GetItem (p->i, (ev->a = 0));
    if (!PyInt_Check(item)) {
      fprintf (stderr, "expecting an integer");
      return -1;
    }
    *i = PyInt_AsLong (item);
  }
  return 0;
}

static int py_end (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  PyObject * item = PyList_GetItem (p->i, ev->a);
  return item != NULL && PyInt_Check(item);
}

static int py_inc (int * i, double * t, Event * ev) {
  PyEvent * p = ev->data;
  PyObject * item = PyList_GetItem (p->i, ++ev->a);
  if (item && PyInt_Check(item))
    *i = PyInt_AsLong (item);
  else
    (*i)++;
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

int py_register_event (PyObject * action, PyObject * i) {
  if (!PyCallable_Check(action)) {
    fprintf (stderr, "parameter must be callable");
    return -1;
  }
  if (!PyList_Check(i) && !PyInt_Check(i)) {
    fprintf (stderr, "parameter must be a list or an int");
    return -1;
  }
  Event ev;
  ev.last = false;
  ev.action = py_action;
  if (PyInt_Check(i)) {
    ev.expr[0] = py_start;
    ev.nexpr = 1;
  }
  else { // list
    ev.expr[0] = py_end; ev.expr[1] = py_start; ev.expr[2] = py_inc;
    ev.nexpr = 3;
  }
  ev.arrayi = NULL;
  ev.arrayt = NULL;
  PyObject * code = PyObject_GetAttrString (action, "__code__");
  ev.file = strdup (PyString_AsString (PyObject_GetAttrString 
				       (code, "co_filename")));
  ev.line = PyInt_AsLong (PyObject_GetAttrString (code, "co_firstlineno"));
  PyObject * name = PyObject_GetAttrString (action, "func_name");
  ev.name = strdup (PyString_AsString (name));
  ev.a = 0;

  // fixme: memory should be freed before exit
  PyEvent * p = ev.data = malloc (sizeof (PyEvent));
  p->i = i;
  Py_INCREF (i);
  p->action = action;
  Py_INCREF (action);

  event_register (ev);
  return 0;
}
