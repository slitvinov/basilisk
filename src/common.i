%{
  typedef int scalar;
  typedef struct {
    scalar x, y;
  } vector;
  typedef struct {
    vector x, y;
  } tensor;

  static void origin (double x, double y) {
    extern double X0, Y0;
    X0 = x; Y0 = y;
  } 
  static void size (double L) {
    extern double L0;
    L0 = L;
  }
  extern void init_solver (void);
  extern void init_grid (int n);
  extern void free_grid (void);
  extern void run (void);

  extern double interpolate (scalar s, double x, double y);
  extern int py_scalar_init (scalar s, PyObject * f);
  extern int py_register_event (PyObject * action, PyObject * i, PyObject * t);
%}

typedef int scalar;

%rename(_vector) vector;
typedef struct {
  scalar x, y;
} vector;

%rename(_tensor) tensor;
typedef struct {
  vector x, y;
} tensor;

extern void origin (double x = 0., double y = 0.);
extern void size (double L = 1.);

extern void init_grid (int n);
extern void free_grid (void);
extern double interpolate (scalar s, double x = 0., double y = 0.);
extern int py_scalar_init (scalar s, PyObject * f);
extern int py_register_event (PyObject * action, PyObject * i, PyObject * t);

%init %{
  init_solver();
%}

%pythoncode %{
class scalar(int):
    def __init__(self,i):
        self = i
    def __setattr__(self, name, value):
        if name == "f":
            py_scalar_init(self,value)
        else:
            self.__dict__[name] = value
    def f(self,x,y=0):
        ret = interpolate(self,x,y)
        if ret > 1e100:
            ret = None
        return ret
    def norm(self):
        return normf(self)
    def stats(self):
        return statsf(self)

class vector:
    def __init__(self,v):
        self.x = scalar(v.x)
        self.y = scalar(v.y)

class tensor:
    def __init__(self,t):
        self.x = vector(t.x)
        self.y = vector(t.y)

def event (action, i = None, t = None):
    py_register_event(action, i, t)
%}
