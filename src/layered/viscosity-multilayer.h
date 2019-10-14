/**
# Viscous friction between layers

Boundary conditions on the top and bottom layers need to be added to close the
system for the viscous stresses. We chose to impose a Neumann condition on the
top boundary i.e.
$$
\partial_z u |_t = \dot{u}_t
$$
and a Navier slip condition on the bottom i.e.
$$
u|_b = u_b + \lambda_b \partial_z u|_b
$$
By default the viscosity is zero and we impose free-slip on the top
boundary and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
$\lambda_b = 0$, $u_b = 0$. */

double nu = 0.;
(const) scalar lambda_b = zeroc, dut = zeroc, u_b = zeroc;

/**
For stability, we discretise the viscous friction term implicitly as
$$
\frac{(hu_l)_{n + 1} - (hu_l)_{\star}}{\Delta t} =
\frac{\nu}{\mathrm{layer}_l}  \left( \frac{u_{l + 1} - u_l}{h_{l + 1 / 2}} -
\frac{u_l - u_{l - 1}}{h_{l - 1 / 2}} \right)_{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Mu}_{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

void vertical_viscosity (Point point, scalar * hl, scalar * sl, double dt)
{
  if (nu == 0.)
    return;
  
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_lu_l = h\mathrm{layer}_lu_l$. */
      
  int l = 0;
  scalar s, h;
  for (s,h in sl,hl)
    rhs[l++] = s[]*h[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{\nu \Delta t}{h_{l - 1 / 2}} \right)_{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{\nu \Delta t}{h_{l + 1 / 2}}
  \right)_{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = \mathrm{layer}_l h_{n + 1} - a_l - c_l
  $$
  */
  
  for (l = 1; l < nl - 1; l++) {
    scalar hm = hl[l-1], h = hl[l], hp = hl[l+1];
    a[l] = - 2.*nu*dt/(hm[] + h[]);
    c[l] = - 2.*nu*dt/(h[] + hp[]);
    b[l] = h[] - a[l] - c[l];
  }
    
  /**
  For the top layer the boundary conditions give the (ghost)
  boundary value
  $$
  u_{\mathrm{nl}} = u_{\mathrm{nl} - 1} + \dot{u}_t h_{\mathrm{nl} - 1},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1} h_{n + 1}
  - a_{\mathrm{nl} - 1}
  $$
  $$
  \mathrm{rhs}_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1}  
  (hu_{\mathrm{nl} - 1})_{\star} + \nu \Delta t \dot{u}_t
  $$
  */

  scalar hm = hl[nl-2]; h = hl[nl-1];
  a[nl-1] = - 2.*nu*dt/(hm[] + h[]);
  b[nl-1] = h[] - a[nl-1];
  rhs[nl-1] += nu*dt*dut[];

  /**
  For the bottom layer, the boundary conditions give the (ghost)
  boundary value $u_{- 1}$
  $$
  u_{- 1} = \frac{2 h_0}{2 \lambda_b + h_0} u_b + \frac{2 \lambda_b - h_0}{2
  \lambda_b + h_0} u_0,
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_0 = \mathrm{layer}_0 h_{n + 1} - c_0 + 
  \frac{2 \nu \Delta t}{2 \lambda_b + h_0}
  $$
  $$
  \mathrm{rhs}_0 = \mathrm{layer}_0  (hu_0)_{\star} + \frac{2 \nu \Delta t}{2
  \lambda_b + h_0} u_b
  $$
  */

#if 0
  scalar hp = hl[1]; h = hl[0];  
  c[0] = - 2.*dt*nu/(h[] + hp[]);
  b[0] = h[] - c[0] + 2.*nu*dt/(2.*lambda_b[] + h[]);
  rhs[0] += 2.*nu*dt/(2.*lambda_b[] + h[])*u_b[];
#else
  scalar h1 = hl[1], h0 = hl[0];
  double den = h[0]*(8.*h[1]*lambda_b[] + sq(h1[])) +
    sq(h0[])*(6.*lambda_b[] + 3.*h[1]) +
    2.*(sq(h[1])*lambda_b[] + cube(h0[]));
  b[0] = h0[] + 2.*dt*nu*(1./(h0[] + h1[]) +
			  (sq(h1[]) + 4.*h[0]*h[1] + 4.*sq(h0[]))/den);
  c[0] = - 2.*dt*nu*(1./(h0[] + h1[]) + sq(h0[])/den);
  rhs[0] += 2.*dt*nu*u_b[]*(sq(h1[]) + 4.*h[0]*h[1] + 3.*sq(h[0]))/den;
#endif
  
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  
  for (l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s = sl[nl-1]; s[] = a[nl-1];
  for (l = nl - 2; l >= 0; l--) {
    s = sl[l];
    s[] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
  }
}
