#define M(i,j) m[(i)*(n + 2) + (j)]
#if 0
  #define NS 3
  #define INDEX(a,k,l) M(i,j).a[k+NS/2][l+NS/2]
#else
  #define NS 1
  #define INDEX(a,k,l) M(i+k,j+l).a[NS/2][NS/2]
#endif

typedef struct {
  double u[NS][NS], v[NS][NS], h[NS][NS], b[NS][NS], ke[NS][NS], psi[NS][NS];
  double un, vn, hn;
} Data;

#define u(k,l)    INDEX(u,k,l)
#define v(k,l)    INDEX(v,k,l)
#define h(k,l)    INDEX(h,k,l)
#define b(k,l)    INDEX(b,k,l)
#define ke(k,l)   INDEX(ke,k,l)
#define psi(k,l)  INDEX(psi,k,l)
#define un(k,l)   M(i+k,j+l).un
#define vn(k,l)   M(i+k,j+l).vn
#define hn(k,l)   M(i+k,j+l).hn

#define foreach(m) for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++)

#define I (i - 1)
#define J (j - 1)

Data * init_grid (int n)
{
  return malloc ((n + 2)*(n + 2)*sizeof (Data));
}

void boundary_h (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++)
	h(k,l) = hn(k,l);

  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    h(0,-1) = h(0,0);
    j = n;
    h(0,1) = h(0,0);
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    h(-1,0) = h(0,0);
    i = n;
    h(1,0) = h(0,0);
  }
}

void boundary_b (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++)
	b(k,l) = M(i+k,j+l).b[NS/2][NS/2];

  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    b(0,-1) = b(0,0);
    j = n;
    b(0,1) = b(0,0);
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    b(-1,0) = b(0,0);
    i = n;
    b(1,0) = b(0,0);
  }
}

void boundary_ke_psi (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++) {
	ke(k,l) = M(i+k,j+l).ke[NS/2][NS/2];
	psi(k,l) = M(i+k,j+l).psi[NS/2][NS/2];
      }

  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    ke(0,-1) = ke(0,0);
    j = n;
    ke(0,1) = ke(0,0);
    psi(0,1) = (v(0,1) - v(-1,1) + u(0,0) - u(0,1))/DX;
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    ke(-1,0) = ke(0,0);
    i = n;
    ke(1,0) = ke(0,0);
    psi(1,0) = (v(1,0) - v(0,0) + u(1,-1) - u(1,0))/DX;
  }
}

void boundary_u (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++) {
	u(k,l) = un(k,l);
	v(k,l) = vn(k,l);
      }

  int i, j;
  for (i = 1; i <= n; i++) {
    /* solid walls */
    j = 1; v(0,0) = 0.;
    j = n; v(0,1) = 0.;
  }
  for (j = 1; j <= n; j++) {
    /* solid walls */
    i = 1; u(0,0) = 0.;
    i = n; u(1,0) = 0.;
  }
}
