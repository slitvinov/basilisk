@include <stdlib.h>
@include <stdio.h>
@include <stddef.h>
@include <stdbool.h>
@include <stdarg.h>
@include <string.h>
@include <float.h>
@include <limits.h>
@include <assert.h>
@include <math.h>
@include <time.h>
@include <sys/time.h>

@define pi 3.14159265358979
@undef HUGE
@define HUGE ((double)1e30)
@define nodata DBL_MAX
@define _NVARMAX 65536
@define is_constant(v) ((v) >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v) - _NVARMAX] : nodata)

@define max(a,b) ((a) > (b) ? (a) : (b))
@define min(a,b) ((a) < (b) ? (a) : (b))
@define sq(x) ((x)*(x))
@define cube(x) ((x)*(x)*(x))
@define sign(x) ((x) > 0 ? 1 : -1)
@define noise() (1. - 2.*rand()/(double)RAND_MAX)
@define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
@define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
@define unmap(x,y)

@define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

@if _OPENMP

@include <omp.h>
@define OMP(x) Pragma(#x)
@define pid() omp_get_thread_num()
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)

@elif _MPI

@include <mpi.h>
@define OMP(x)

@if TRACE
   FILE * tracefp = NULL;
   static double tracestart = 0.;
   static void trace (const char * format, ...) {
     va_list ap;
     va_start (ap, format);
     vfprintf (tracefp, format, ap);
     va_end (ap);
   }
@else
@  define trace(...)
@endif

static bool in_prof = false;
@def prof_start(name)
  assert (!in_prof); in_prof = true;
  double prof_start = MPI_Wtime();
  trace ("10 %.6f State Thread%d \"%s\"\n",
	 prof_start - tracestart, pid(), name)
@
@def prof_stop()
  assert (in_prof); in_prof = false;
  double _prof = MPI_Wtime();
  mpi_time += _prof - prof_start;
  trace ("10 %.6f State Thread%d compute\n", _prof - tracestart, pid())
@
@def trace_event(name)
  trace ("11 %.6f Event Thread%d \"%s\"\n",
	 MPI_Wtime() - tracestart, pid(), name)
@

@if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)
@else
@def mpi_all_reduce(v,type,op) {
  prof_start ("mpi_all_reduce");
  union { int a; float b; double c;} global;
  MPI_Allreduce (&(v), &global, 1, type, op, MPI_COMM_WORLD);
  memcpy (&(v), &global, sizeof (v));
  prof_stop();
}
@
@def mpi_all_reduce_double(v,op) {
  prof_start ("mpi_all_reduce");
  double global, tmp = v;
  MPI_Allreduce (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);
  v = global;
  prof_stop();
}
@

@endif

static int mpi_rank, mpi_npe;
@define pid() mpi_rank
@define npe() mpi_npe

static void finalize (void)
{
@if TRACE
  double t = MPI_Wtime() - tracestart;
  fprintf (tracefp, "8 %.8f Thread%d Thread\n", t, pid());
  if (pid() == 0)
    fprintf (tracefp, "8 %.8f C_Prog Prog\n", t);
  fclose (tracefp);
@endif
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
    char name[80];
    if (mpi_rank > 0) {
      sprintf (name, "out-%d", mpi_rank);
      stdout = freopen (name, "w", stdout);
      sprintf (name, "log-%d", mpi_rank);
      stderr = freopen (name, "w", stderr);
    }
@if TRACE
    sprintf (name, "trace-%d", mpi_rank);
    tracefp = fopen (name, "w");
    if (pid() == 0) {
      #include "paje.h"
    }
    fprintf (tracefp, "7 0 Thread%d Thread C_Prog \"Thread %d\"\n",
	     mpi_rank, mpi_rank);
    fprintf (tracefp, "10 0 State Thread%d compute\n", mpi_rank);
    MPI_Barrier (MPI_COMM_WORLD);
    tracestart = MPI_Wtime();
    double global;
    MPI_Allreduce (&tracestart, &global, 1,
		   MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    tracestart = global;
@endif
  }
}

@else // not MPI, not OpemMP

@define OMP(x)
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)

@endif

// fixme: _OMPSTART and _OMPEND are only used for working around the
// lack of min|max reduction operations in OpenMP < 3.1
@define OMP_PARALLEL()     OMP(omp parallel) { _OMPSTART
@define OMP_END_PARALLEL() _OMPEND }
@define _OMPSTART
@define _OMPEND

@define NOT_UNUSED(x) (void)(x)

@define VARIABLES      _CATCH;
@define val(a,k,l)     data(k,l)[a]

/* undefined value */
/* Initialises unused memory with "signaling NaNs".  
 * This is probably not very portable, tested with
 * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
 * This blog was useful:
 *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
 */
@if _GNU_SOURCE
double undefined;
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  feenableexcept (FE_DIVBYZERO|FE_INVALID);
}
@else
@  define undefined DBL_MAX
@endif

// the grid
void * grid = NULL;
// coordinates of the lower-left corner of the box
static double X0 = 0., Y0 = 0.;
// size of the box
static double L0 = 1.;
// number of grid points
static int N = 64;

typedef int scalar;

typedef struct {
  scalar x, y;
} vector;

typedef struct {
  vector x, y;
} tensor;

#define norm(v) (sqrt(sq(v.x[]) + sq(v.y[])))

void origin (double x, double y) {
  X0 = x; Y0 = y;
}

void size (double L) {
  L0 = L;
}

static double nothing() { return 0.; }
void * none = nothing;

// boundary conditions for each direction/variable

enum { right, left, top, bottom, nboundary };
// ghost cell coordinates for each direction
int _ig[nboundary] = {1,-1,0,0}, 
    _jg[nboundary] = {0,0,1,-1};

@define dirichlet(x)            (2.*(x) - val(_s,0,0))
@define dirichlet_homogeneous() (- val(_s,0,0))
@define neumann(x)              (Delta*(x) + val(_s,0,0))
@define neumann_homogeneous()   (val(_s,0,0))

double  * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#include "grid/boundaries.h"

// attributes for each scalar

@include "_attributes.h"

attribute {
  double (* boundary[nboundary])             (Point, scalar);
  double (* boundary_homogeneous[nboundary]) (Point, scalar);
  double (* gradient)                        (double, double, double);
  char * name;
  struct { int x, y; } d; // staggering
  vector v;
  bool   face, normal;
};

// lists

int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = realloc (list, sizeof (scalar)*(len + 2));
  list[len] = s;
  list[len + 1] = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t == s)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    for (scalar s1 in l)
      if (s1 == s)
	return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    for (scalar s in l)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  for (scalar s in l2)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  for (scalar s in l)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", s.name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

/**
Given a list of scalars this functions splits it into two lists: a
list of scalars defined on the faces in direction d and the
rest. */

void list_split (scalar * list, int d,
		 scalar ** faces, scalar ** rest)
{
  *faces = *rest = NULL;
  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[d]) {
      if (s.face && d % 2 && (&s.d.x)[component])
	*faces = list_add (*faces, s);
      else
	*rest = list_add (*rest, s);
    }
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = realloc (list, sizeof (vector)*(len + 2));
  list[len] = v;
  list[len + 1] = (vector){-1,-1};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  for (vector w in list)
    if (w.x == v.x && w.y == v.y)
      return list;
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    for (vector v in l)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (*s >= 0) {
    vector v;
    v.x = *s++;
    assert (*s >= 0);
    v.y = *s++;
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  for (tensor t in list) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = realloc (list, sizeof (tensor)*(len + 2));
  list[len] = t;
  list[len + 1] = (tensor){{-1,-1},{-1,-1}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x >= 0) {
    tensor t;
    t.x = *v++;
    assert (v->x >= 0);
    t.y = *v++;
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL; // all the fields

// basic methods

scalar (* init_scalar)      (scalar, const char *);
vector (* init_vector)      (vector, const char *);
tensor (* init_tensor)      (tensor, const char *);
vector (* init_face_vector) (vector, const char *);

// events 

typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
};

static Event * Events = NULL; // all events

double tnext = 0; // time of next event
void init_events (void);
void event_register (Event event);

void init_solver (void);

// timers

@if _MPI
static double mpi_time = 0.;
@endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
@if _MPI
  t.tm = mpi_time;
@endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) + 
	  (tvend.tv_usec - t.tv.tv_usec)/1e6);
}
	  
// Constant fields

const face vector zerof[] = {0.,0.};
const face vector unityf[] = {1.,1.};
const scalar unity[] = 1.;

// Metric

(const) face vector fm = unityf;
(const) scalar cm = unity;

// Arrays

typedef struct {
  void * p;
  size_t size, max, len;
} Array;

Array * array_new (size_t size)
{
  Array * a = malloc (sizeof(Array));
  a->p = NULL;
  a->size = size;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  if (a->max > 0)
    free (a->p);
  free (a);
}

void array_append (Array * a, void * elem)
{
  if (a->len == a->max) {
    a->max += 128;
    a->p = realloc (a->p, a->max*a->size);
  }
  memcpy (a->p + a->len++*a->size, elem, a->size);
}

void array_swap (Array * a, int i, int j)
{
  char buf[a->size];
  memcpy (buf, a->p + i*a->size, a->size);
  memcpy (a->p + i*a->size, a->p + j*a->size, a->size);
  memcpy (a->p + j*a->size, buf, a->size);
}

void array_reverse (Array * a)
{
  for (int i = 0; i < a->len/2; i++)
    array_swap (a, i, a->len - 1 - i);
}
