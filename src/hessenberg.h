#define sign(a,b) ((b) > 0. ? ((a) > 0. ? (a) : -(a)) : ((a) > 0. ? -(a) : (a)))

static inline void givens (double x, double y, double * c, double * s)
{
#if 0
  if (x == 0. && y == 0.)
    *c = 1., *s = 0.;
  else if (fabs(y) > fabs(x)) {
    double t = x/y;
    x = sqrt(1. + t*t);
    *s = - sign(1./x, y);
    *c = t*(*s);
  }
  else {
    double t = y/x;
    y = sqrt(1. + t*t);
    *c = sign(1./y, x);
    *s = - t*(*c);
  }
#else
  double t = sqrt (sq(x) + sq(y));
  *c = x/t, *s = -y/t;
#endif
}

/** Henry 1994 algorithm */
void solve_hessenberg (double * H, double * x, int n)
{
  double v[n], c[n], s[n];
  for (int i = 0; i < n; i++)
    v[i] = H[n*(i + 1) - 1];
  for (int k = n - 1; k >= 1; k--) {
    double a = H[k*n + k - 1];
    givens (v[k], a, &c[k], &s[k]);
    x[k] /= c[k]*v[k] - s[k]*a;
    double ykck = x[k]*c[k], yksk = x[k]*s[k];
    for (int l = 0; l <= k - 2; l++) {
      a = H[l*n + k - 1];
      x[l] -= ykck*v[l] - yksk*a;
      v[l] = c[k]*a + s[k]*v[l];
    }
    a = H[(k - 1)*n + k - 1];
    x[k-1] -= ykck*v[k-1] - yksk*a;
    v[k-1] = c[k]*a + s[k]*v[k-1];
  }
  double tau1 = x[0]/v[0];
  for (int k = 1; k < n; k++) {
    double tau2 = x[k];
    x[k-1] = c[k]*tau1 - s[k]*tau2;
    tau1 = c[k]*tau2 + s[k]*tau1;
  }
  x[n-1] = tau1;
}
