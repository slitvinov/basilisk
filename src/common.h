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
@include <sys/resource.h>

@define pi 3.14159265358979
@undef HUGE
@define HUGE ((double)1e30)
@define nodata HUGE
@define _NVARMAX 65536
@define is_constant(v) ((v).i >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

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

@define systderr  stderr
@define systdout  stdout

@if _MPI
static FILE * qstderr (void);
static FILE * qstdout (void);
FILE * ferr, * fout;
@else
@ define qstderr() stderr
@ define qstdout() stdout
@ define ferr      stderr
@ define fout      stdout
@endif

// Memory tracing

@if MTRACE

struct {
  FILE * fp;                     // trace file
  size_t total, max;             // current and maximum allocated memory
  size_t overhead, maxoverhead;  // current and maximum profiling overhead
  size_t nr;                     // current number of records
  size_t startrss, maxrss;       // starting and maximum system ressource usage
  char * fname;                  // trace file name
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

@define sysmalloc malloc
@define syscalloc calloc
@define sysrealloc realloc
@define sysfree free
@define systrdup strdup

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
	     c, f->id, pmtrace.nr, pmtrace.total, f->total);
@if (_GNU_SOURCE || _DARWIN_C_SOURCE)
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
@endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
			    const char * func, const char * file, int line,
			    char c)
{
  assert (d != NULL);
  d->id = pmfunc_index(func, file, line);
  d->size = size;
  pmfunc * f = &pmfuncs[d->id - 1];
  f->total += size;
  if (f->total > f->max)
    f->max = f->total;
  pmtrace.total += size;
  pmtrace.overhead += sizeof(pmdata);
  if (pmtrace.total > pmtrace.max) {
    pmtrace.max = pmtrace.total;
    pmtrace.maxoverhead = pmtrace.overhead;
  }
  pmfunc_trace (f, c);
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", stderr);
    if (d->size == 0)
      fputs (", possible double free()", stderr);
    else
      fputs (", not traced?", stderr);
    fputs (", aborting...\n", stderr);
    abort();
    return ptr;
  }
  else {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (stderr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
	       f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (stderr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
	       pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
    return d;
  }
}

static void * pmalloc (size_t size,
		       const char * func, const char * file, int line)
{
  return pmfunc_alloc (sysmalloc (sizeof(pmdata) + size),
		       size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
		       const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
			const char * func, const char * file, int line)
{
  return pmfunc_alloc (sysrealloc (pmfunc_free(ptr, '<'),
				   sizeof(pmdata) + size),
		       size, func, file, line, '>');
}

static void pfree (void * ptr,
		   const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
		       const char * func, const char * file, int line)
{
  char * d = pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

@if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
@endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
@if MTRACE < 3
  fprintf (stderr,
	   "*** MTRACE: max resident  set size: %10ld bytes\n"
	   "*** MTRACE: max traced memory size: %10ld bytes"
	   " (tracing overhead %.1g%%)\n"
	   "%10s    %20s   %s\n",
	   pmtrace.maxrss*1024,
	   pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
	   "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (stderr, "%10ld    %20s   %s:%d\n",
	     p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
	     "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
	     "total(\"%s\") w l t 'total'",
	     fname,
	     pmtrace.startrss*1024.,
	     pmtrace.startrss*1024.,
	     fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
	       ",func(\"%s\",%d) w l t '%s'",
	       fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (stderr,
	     "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
	     fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
@endif // MTRACE < 3

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (stderr, "%s:%d: error: %ld bytes leaked here\n",
	       p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
@if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", stderr);
@endif
    pmfuncs_free();
  }
}

@else // !MTRACE
@ define pmalloc(s,func,file,line)    malloc(s)
@ define pcalloc(n,s,func,file,line)  calloc(n,s)
@ define prealloc(p,s,func,file,line) realloc(p,s)
@ define pfree(p,func,file,line)      free(p)
@ define pstrdup(s,func,file,line)    strdup(s)
@endif // !MTRACE

// Arrays

typedef struct {
  char * p;
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

// Function tracing

@if TRACE
@include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func     = {
  {NULL, sizeof(char *), 0, 0}, {NULL, sizeof(int), 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, sizeof(char *), 0, 0}, {NULL, sizeof(int), 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len; i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = strdup (func);
  array_append (a, &s);
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value);
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len--;
  int value = t->stack.len > 0 ? ((int *)t->stack.p)[t->stack.len - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len + 1];
    char * names[t->index.len + 1], ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len; i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len; i++, func++)
    free (*func);
  free (t->index.p);
  free (t->stack.p);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}
#if 0
#define TRACE_TYPE(func) (strncmp (func, "mpi_", 4) ?		\
			  &trace_func : &trace_func)
#else
#define TRACE_TYPE(func) &trace_func
#endif
@  define trace(func, file, line)     trace_push (TRACE_TYPE(func), func)
@  define end_trace(func, file, line) trace_pop (TRACE_TYPE(func), func)
@else // !TRACE
@  define trace(...)
@  define end_trace(...)
@endif

// OpenMP / MPI
  
@if _OPENMP

@include <omp.h>
@define OMP(x) Pragma(#x)
@define tid() omp_get_thread_num()
@define pid() 0
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)

@elif _MPI

@include <mpi.h>
@define OMP(x)

static bool in_prof = false;
static double prof_start, _prof;
@def prof_start(name)
  assert (!in_prof); in_prof = true;
  prof_start = MPI_Wtime();
@
@def prof_stop()
  assert (in_prof); in_prof = false;
  _prof = MPI_Wtime();
  mpi_time += _prof - prof_start;
@

@if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)
@else
trace
static int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
			    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);
}
@def mpi_all_reduce(v,type,op) {
  prof_start ("mpi_all_reduce");
  union { int a; float b; double c;} global;
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);
  memcpy (&(v), &global, sizeof (v));
  prof_stop();
}
@
@def mpi_all_reduce_double(v,op) {
  prof_start ("mpi_all_reduce");
  double global, tmp = v;
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);
  v = global;
  prof_stop();
}
@

@endif

static int mpi_rank, mpi_npe;
@define tid() mpi_rank
@define pid() mpi_rank
@define npe() mpi_npe

@define QFILE FILE // a dirty trick to avoid qcc 'static FILE *' rule

static FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

static FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
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
    srand (mpi_rank + 1);
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = stderr;
      fout = stdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
@if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
@endif
  }
}

@else // not MPI, not OpenMP

@define OMP(x)
@define tid() 0
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_double(v,op)

@endif // not MPI, not OpenMP

void init_solver()
{
@if _MPI
  mpi_init();
@elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
@endif
}

// fixme: _OMPSTART and _OMPEND are only used for working around the
// lack of min|max reduction operations in OpenMP < 3.1
@define OMP_PARALLEL()     OMP(omp parallel) { _OMPSTART
@define OMP_END_PARALLEL() _OMPEND }
@define _OMPSTART
@define _OMPEND

@define NOT_UNUSED(x) (void)(x)

@define VARIABLES      _CATCH;
@define val(a,k,l,m)   data(k,l,m)[a.i]

double _val_higher_dimension = 0.;
@define _val_higher_dimension(x,a,b,c) _val_higher_dimension

/* undefined value */
/* Initialises unused memory with "signaling NaNs".  
 * This is probably not very portable, tested with
 * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
 * This blog was useful:
 *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
 */
@if (_GNU_SOURCE || __APPLE__)
double undefined;
@ if __APPLE__
@   include <stdint.h>
@   include "fp_osx.h"
@ endif
@  define enable_fpe(flags)  feenableexcept (flags)
@  define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
@else
@  define undefined DBL_MAX
@  define enable_fpe(flags)
@  define disable_fpe(flags)
static void set_fpe (void) {}
@endif

// the grid
void * grid = NULL;
// coordinates of the lower-left corner of the box
double X0 = 0., Y0 = 0., Z0 = 0.;
// size of the box
double L0 = 1.;
// number of grid points
int N = 64;

typedef struct { int i; } scalar;

typedef struct {
  scalar x;
#if dimension > 1
  scalar y;
#endif
#if dimension > 2
  scalar z;
#endif
} vector;

typedef struct {
  vector x;
#if dimension > 1
  vector y;
#endif
#if dimension > 2
  vector z;
#endif
} tensor;

typedef struct {
  double x, y, z;
} coord;

#if dimension == 1
# define norm(v) fabs(v.x[])
# define dv() (Delta*cm[])
#elif dimension == 2
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[])))
# define dv() (sq(Delta)*cm[])
#else // dimension == 3
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[]) + sq(v.z[])))
# define dv() (cube(Delta)*cm[])
#endif

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }

// boundary conditions for each direction/variable

#if dimension == 1
  enum { right, left };
#elif dimension == 2
  enum { right, left, top, bottom };
#else
  enum { right, left, top, bottom, front, back };
#endif
int nboundary = 2*dimension;

#define none -1

@define dirichlet(x)            (2.*(x) - val(_s,0,0,0))
@define dirichlet_homogeneous() (- val(_s,0,0,0))
@define neumann(x)              (Delta*(x) + val(_s,0,0,0))
@define neumann_homogeneous()   (val(_s,0,0,0))

double  * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#include "grid/boundaries.h"

// attributes for each scalar

@include "_attributes.h"

attribute {
  double (** boundary)             (Point, Point, scalar);
  double (** boundary_homogeneous) (Point, Point, scalar);
  double (* gradient)              (double, double, double);
  void   (* delete)                (scalar);
  char * name;
  struct {
    int x;
#if dimension > 1
    int y;
#endif
#if dimension > 2
    int z;
#endif
  } d; // staggering
  vector v;
  bool   face, normal;
}

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
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    for (scalar s1 in l)
      if (s1.i == s.i)
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
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  for (vector w in list) {
    bool id = true;
    foreach_dimension()
      if (w.x.i != v.x.i)
	id = false;
    if (id)
      return list;
  }
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
  while (s->i >= 0) {
    vector v;
    foreach_dimension() {
      assert (s->i >= 0);
      v.x = *s++;
    }
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
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    foreach_dimension() {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
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

#define vector(x) (*((vector *)&(x)))

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
void _init_solver (void);

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

const face vector zerof[] = {0.,0.,0.};
const face vector unityf[] = {1.,1.,1.};
const scalar unity[] = 1.;
const scalar zeroc[] = 0.;

// Metric

(const) face vector fm = unityf;
(const) scalar cm = unity;

// Pipes

static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = realloc (qpopen_pipes, sizeof(FILE *)*(n + 2));
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  free (qpopen_pipes);
  qpopen_pipes = NULL;
}

#define popen  qpopen
#define pclose qpclose
