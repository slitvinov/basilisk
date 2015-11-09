#include "utils.h"

#define PARABOLA_FIT_CENTER_WEIGHT .1
#define PARABOLA_SIMPLER 0

typedef struct {
  coord o;
#if dimension == 2 /* y = a[0]*x^2 + a[1]*x + a[2] */
  coord m;
  double ** M;
  double rhs[3], a[3];
#else /* 3D */
# if PARABOLA_SIMPLER /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y */
  GtsMatrix * M;
  GtsVector rhs, a;
# else /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y + a[3]*x + a[4]*y + a[5] */
  gdouble ** M, rhs[6], a[6];
# endif
  GtsVector t[3];
#endif /* 3D */
} ParabolaFit;

static void normalize (coord * n)
{
  double norm = 0.;
  foreach_dimension()
    norm += sq(n->x);
  norm = sqrt(norm);
  foreach_dimension()
    n->x /= norm;
}

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  foreach_dimension()
    p->o.x = o.x;
#if dimension == 2
  foreach_dimension()
    p->m.x = m.x;
  normalize (&p->m);
  p->M = matrix_new (3, 3, sizeof(double));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
#else /* 3D */
  gdouble max;
  GtsVector nx = {0., 0., 0.}, ny, nz;
  guint d = 0;

  nz[0] = m->x; nz[1] = m->y; nz[2] = m->z;
  gts_vector_normalize (nz);
  max = nz[0]*nz[0];
  /* build a vector orthogonal to nz */
  if (nz[1]*nz[1] > max) { max = nz[1]*nz[1]; d = 1; }
  if (nz[2]*nz[2] > max) d = 2;
  switch (d) {
  case 0: nx[0] = - nz[2]/nz[0]; nx[2] = 1.0; break;
  case 1: nx[1] = - nz[2]/nz[1]; nx[2] = 1.0; break;
  case 2: nx[2] = - nz[0]/nz[2]; nx[0] = 1.0; break;
  }
  gts_vector_normalize (nx);

  /* build a second vector orthogonal to nx and nz */
  gts_vector_cross (ny, nz, nx);

  /* transformation matrix from (i,j,k) to (nx, ny, nz) */
  p->t[0][0] = nx[0]; p->t[0][1] = nx[1]; p->t[0][2] = nx[2];
  p->t[1][0] = ny[0]; p->t[1][1] = ny[1]; p->t[1][2] = ny[2];
  p->t[2][0] = nz[0]; p->t[2][1] = nz[1]; p->t[2][2] = nz[2];

# if PARABOLA_SIMPLER
  p->M = gts_matrix_zero (NULL);
  p->rhs[0] = p->rhs[1] = p->rhs[2] = 0.;
# else
  p->M = gfs_matrix_new (6, 6, sizeof (gdouble));
  p->rhs[0] = p->rhs[1] = p->rhs[2] = p->rhs[3] = p->rhs[4] = p->rhs[5] = 0.;
# endif
#endif /* 3D */
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{
#if dimension == 2
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#else /* 3D */
  gdouble x1 = m[0] - p->o[0];
  gdouble y1 = m[1] - p->o[1];
  gdouble z1 = m[2] - p->o[2];
  gdouble x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  gdouble y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  gdouble z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
  gdouble x2 = x*x, x3 = x2*x, x4 = x3*x;
  gdouble y2 = y*y, y3 = y2*y, y4 = y3*y;
# if PARABOLA_SIMPLER
  p->M[0][0] += w*x4;
  p->M[1][0] += w*x2*y2; p->M[1][1] += w*y4;
  p->M[2][0] += w*x3*y;  p->M[2][1] += w*x*y3;
  p->rhs[0] += w*z*x2;   p->rhs[1] += w*z*y2;  p->rhs[2] += w*z*x*y;
# else
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2; 
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;
# endif
#endif /* 3D */
}

static void parabola_fit_solve (ParabolaFit * p)
{
#if dimension == 2
  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  if (matrix_inverse (p->M, 3, 1e-6)) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = 0.;
#else /* 3D */
# if PARABOLA_SIMPLER
  p->M[0][1] = p->M[1][0]; p->M[0][2] = p->M[2][0];
  p->M[1][2] = p->M[2][1]; p->M[2][2] = p->M[1][0];
  GtsMatrix * M = gts_matrix3_inverse ((GtsMatrix *) p->M);
  if (M) {
    p->a[0] = M[0][0]*p->rhs[0] + M[0][1]*p->rhs[1] + M[0][2]*p->rhs[2];
    p->a[1] = M[1][0]*p->rhs[0] + M[1][1]*p->rhs[1] + M[1][2]*p->rhs[2];
    p->a[2] = M[2][0]*p->rhs[0] + M[2][1]*p->rhs[1] + M[2][2]*p->rhs[2];
    gts_matrix_destroy (M);
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = p->a[2] = 0.;
# else
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  guint i, j;
  for (i = 1; i < 6; i++)
    for (j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  if (gfs_matrix_inverse (p->M, 6, 1e-10)) {
    for (i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (j = 0; j < 6; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  }
  else { /* this may be a degenerate/isolated interface fragment */
    g_warning ("singular matrix");
    p->a[0] = p->a[1] = p->a[2] = 0.;
  }
# endif
#endif /* 3D */  
  matrix_free (p->M);
}

static double parabola_fit_curvature (ParabolaFit * p,
				      double kappamax, double * kmax)
{
  double kappa;
#if dimension == 2
  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#else /* 3D */
  gdouble hxx = 2.*p->a[0];
  gdouble hyy = 2.*p->a[1];
  gdouble hxy = p->a[2];
  gdouble hx, hy;
# if PARABOLA_SIMPLER
  hx = hy = 0.;
# else
  hx = p->a[3];
  hy = p->a[4];
# endif
  gdouble dnm = 1. + hx*hx + hy*hy;
  kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/sqrt (dnm*dnm*dnm);
  if (kmax) {
    gdouble kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    gdouble a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }
#endif /* 3D */
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}

#if AXI
static void parabola_fit_axi_curvature (const ParabolaFit * p,
					double r, double h,
					double * kappa, double * kmax)
{
  double nr = (p->m.x*p->a[1] + p->m.y)/sqrt (1. + sq(p->a[1]));
  /* limit the minimum radius to half the grid size */
  double kaxi = nr/max(r, h/2.);
  *kappa += kaxi;
  if (kmax)
    *kmax = max (*kmax, fabs (kaxi));
}
#endif /* 2D */
