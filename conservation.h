// generic solver for system of conservation laws
#include "utils.h"

typedef struct {
  double l, r;
} state; // left-right Riemann states

// primary variables (i.e. will be coarsened/refined)
scalar * primary = NULL;

// extra storage for predictor/corrector
// fixme: not needed for first-order scheme
scalar * conserved1 = NULL;

vector * slopes = NULL;
scalar * tendencies = NULL;

// Default user-provided parameters
// conserved quantities
scalar * conserved = NULL;
// gradient
double (* gradient) (double, double, double) = minmod2;
// user-provided functions
void   parameters   (void);
void   init         (void);
double riemann      (state * s, double delta, double * flux, double dtmax);

int list_len (scalar * list)
{
  int len = 0;
  for (scalar s in list)
    len++;
  return len;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * list = malloc (sizeof(scalar)*(list_len(l1) + list_len(l2) + 1));
  int i = 0;
  for (scalar s in l1)
    list[i++] = s;
  for (scalar s in l2)
    list[i++] = s;
  list[i] = -1;
  return list;
}

static double flux (double dtmax)
{
  int len = list_len (conserved);
  state  c[len];
  double f[len];

  foreach_face() {
    double dx = delta/2.;
    int i = 0;
    scalar s;
    vector g;
    for (s,g in conserved,slopes) {
      c[i].l = s[] - dx*g.x[];
      c[i++].r = s[-1,0] + dx*g.x[-1,0];
    }
    // Riemann solver
    dtmax = riemann (c, delta, f, dtmax);
    // update tendencies
    i = 0;
    for (scalar ds in tendencies) {
      ds[] += f[i];
      ds[-1,0] -= f[i++];
    }
  }

#if QUADTREE
  // propagate updates from fine to coarse
  foreach_halo()
    for (scalar ds in tendencies) {
      coarse(ds,0,0) += ds[]/2.;
      ds[] = 0.;
    }
#endif

  return dtmax;
}

static void update (scalar * list1, scalar * list2, double dt)
{
  if (list1 != list2)
    trash (list1);
  foreach() {
    scalar s1, s2, ds;
    for (s1,s2,ds in list1,list2,tendencies) {
      s1[] = s2[] + dt*ds[]/delta;
      ds[] = 0.;
    }
  }
  boundary (list1);
}

double dt = 0.;

void run()
{
  parameters();
  init_grid(N);

#if QUADTREE
  // we need the tendencies to be reinitialised during refinement
  for (scalar ds in tendencies)
    ds.refine = refine_reset;
#endif

  // limiting
  for (scalar s in conserved)
    s.gradient = gradient;
  scalar * list = (scalar *) slopes;
  for (scalar s in list)
    s.gradient = zero;

  // default values
  primary = list_concat (conserved, tendencies);
  foreach()
    for (scalar s in primary)
      s[] = 0.;
  boundary (primary);

  // user-defined initial conditions
  init();
  boundary (primary);

  // clone temporary storage
  clone (conserved, conserved1);

  // main loop
  timer start = timer_start();
  double t = 0.;
  int i = 0, tnc = 0;
  while (events (i, t)) {
    gradients (conserved, slopes);
    dt = dtnext (t, flux (DT));

    if (gradient == zero)
      /* 1st-order time-integration for 1st-order spatial
	 discretisation */
      update (conserved, conserved, dt);
    else {
      /* 2nd-order time-integration */
      /* predictor */
      update (conserved1, conserved, dt/2.);
      swap (scalar *, conserved, conserved1);
      
      /* corrector */
      gradients (conserved, slopes);
      flux (dt);

      update (conserved, conserved1, dt);
    }

    foreach (reduction(+:tnc))
      tnc++;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free (primary);
  free_grid ();
}
