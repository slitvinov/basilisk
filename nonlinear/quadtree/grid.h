#include "quadtree.c"

#define M(k,l) m[quad_neighbor_finest(p, k, l)]

typedef struct {
  double u, v, h, b, ke, psi;
  double un, vn, hn;
} Data;

#define u(k,l)    M(k,l).u
#define v(k,l)    M(k,l).v
#define h(k,l)    M(k,l).h
#define b(k,l)    M(k,l).b
#define ke(k,l)   M(k,l).ke
#define psi(k,l)  M(k,l).psi
#define un(k,l)   M(k,l).un
#define vn(k,l)   M(k,l).vn
#define hn(k,l)   M(k,l).hn

#define foreach(m) for (int p = size(quad_r - 1), stop = size(quad_r); p < stop; p++)

#define I quad_x(p)
#define J quad_y(p)

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

  /* periodic by construction */
}

void boundary_u (Data * m, int n)
{
  /* update stencils */
  foreach (m) {
    u(0,0) = un(0,0);
    v(0,0) = vn(0,0);
  }

  /* periodic by construction */
}

void boundary_b (Data * m, int n)
{
  /* update stencils */
  
  /* periodic by construction */
}

void boundary_ke_psi (Data * m, int n)
{
  /* update stencils */

  /* periodic by construction */
}
