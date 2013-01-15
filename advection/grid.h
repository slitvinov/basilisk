#include <stdio.h>
#include "mgrid.c"

#define M(i,j) m[(i)*(n + 2) + (j)]

typedef struct {
  double u, v, h, b, ke, psi;
  double un, vn, hn;
} Data;

#define u(k,l)    M(i+k,j+l).u
#define v(k,l)    M(i+k,j+l).v
#define h(k,l)    M(i+k,j+l).h
#define b(k,l)    M(i+k,j+l).b
#define ke(k,l)   M(i+k,j+l).ke
#define psi(k,l)  M(i+k,j+l).psi
#define un(k,l)   M(i+k,j+l).un
#define vn(k,l)   M(i+k,j+l).vn
#define hn(k,l)   M(i+k,j+l).hn

#define I (i - 1)
#define J (j - 1)

Data * init_grid (int n)
{
  int r = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    r++;
  }
  return quadtree (r, sizeof (Data));
}

void boundary_h (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    h(0,0) = hn(0,0);

  for (int i = 1; i <= n; i++) {
    /* periodic */
    M(i,0).h = M(i,n).h;
    M(0,i).h = M(n,i).h;
    M(i,n+1).h = M(i,1).h;
    M(n+1,i).h = M(1,i).h;
  }
}

void boundary_b (Data * m, int n)
{
  /* update stencils */

  for (int i = 1; i <= n; i++) {
    /* periodic */
    M(i,0).b = M(i,n).b;
    M(0,i).b = M(n,i).b;
    M(i,n+1).b = M(i,1).b;
    M(n+1,i).b = M(1,i).b;
  }
}

void boundary_ke_psi (Data * m, int n)
{
  /* update stencils */

  for (int i = 1; i <= n; i++) {
    /* periodic */
    M(i,0).ke = M(i,n).ke;
    M(0,i).ke = M(n,i).ke;
    M(i,n+1).ke = M(i,1).ke;
    M(n+1,i).ke = M(1,i).ke;    
    /* periodic */
    M(i,0).psi = M(i,n).psi;
    M(0,i).psi = M(n,i).psi;
    M(i,n+1).psi = M(i,1).psi;
    M(n+1,i).psi = M(1,i).psi;    
  }
}

void boundary_u (Data * m, int n)
{
  /* update stencils */
  foreach (m) {
    u(0,0) = un(0,0);
    v(0,0) = vn(0,0);
  }

  for (int i = 1; i <= n; i++) {
    /* periodic */
    M(i,0).u = M(i,n).u;
    M(0,i).u = M(n,i).u;
    M(i,n+1).u = M(i,1).u;
    M(n+1,i).u = M(1,i).u;    
    /* periodic */
    M(i,0).v = M(i,n).v;
    M(0,i).v = M(n,i).v;
    M(i,n+1).v = M(i,1).v;
    M(n+1,i).v = M(1,i).v;
  }
}
