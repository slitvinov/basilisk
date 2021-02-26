/**
# Non-hydrostatic extension of the multilayer solver

This adds the non-hydrostatic terms of the [vertically-Lagrangian
multilayer solver for free-surface flows](hydro.h) described in
[Popinet, 2020](/Bibliography#popinet2020). The corresponding system
of equations is
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta) 
  {\color{blue} - \mathbf{{\nabla}} (h \phi)_k + \left[ \phi 
      \mathbf{{\nabla}} z \right]_k},\\
  {\color{blue} \partial_t (hw)_k + \mathbf{{\nabla}} \cdot 
  \left( hw \mathbf{u} \right)_k} & {\color{blue} = - [\phi]_k,}\\
  {\color{blue} \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k + 
  \left[ w - \mathbf{u} \cdot \mathbf{{\nabla}} z \right]_k} & 
  {\color{blue} = 0},
\end{aligned}
$$
where the terms in blue are non-hydrostatic.

The additional $w_k$ and $\phi_k$ fields are defined. The convergence
statistics of the multigrid solver are stored in *mgp*.

Wave breaking is parameterised usng the *breaking* parameter, which is
turned off by default (see Section 3.6.4 in [Popinet,
2020](/Bibliography#popinet2020)). */

#define NH 1
#include "implicit.h"

scalar w, q;
mgstats mgp;
double breaking = HUGE;

/**
## Setup

The $w_k$ and $\phi_k$ scalar fields are allocated and the $w_k$ are
added to the list of advected tracers. */

event defaults (i = 0)
{
  hydrostatic = false;
  if (CFL_H == nodata)
    CFL_H = 1.; // fixme: HUGE as set in implicit.h
  mgp.nrelax = 4;
  
  assert (nl > 0);
  w = new scalar[nl];
  q = new scalar[nl];
  reset ({w, q}, 0.);

  if (!linearised)
    tracers = list_append (tracers, w);
}

/**
## Viscous term

Vertical diffusion is added to the vertical component of velocity
$w$. */

event viscous_term (i++)
{
  if (nu > 0.) {
    // fixme: ugly hack
    scalar lb = lambda_b, d = dut, u = u_b;
    lambda_b = dut = u_b = zeroc;
    foreach()
      vertical_viscosity (point, h, w, dt);
    boundary ({w});
    lambda_b = lb; dut = d; u_b = u;
  }
}

/**
## Assembly of the Hessenberg matrix

For the Keller box scheme, the linear system of equations verified by
the non-hydrostatic pressure $\q$ is expressed as an [Hessenberg
matrix](https://en.wikipedia.org/wiki/Hessenberg_matrix) for each column.

The Hessenberg matrix $\mathbf{H}$ for a column at a particular *point* is
stored in a one-dimensional array with `nl*nl` elements. It encodes
the coefficients of the left-hand-side of the Poisson equation as
$$
\begin{aligned}
  (\mathbf{H}\mathbf{\q} - \mathbf{d})_l & =
  - \text{rhs}_l +
  h_l \partial_x \partial_x (h_l \q_{l - 1 / 2}) + 
  h_l \partial_x \partial_x (h_l \q_{l + 1 / 2}) +\\
  & 4 (\q_{l + 1 / 2} - \q_{l - 1 / 2}) + 8 h_l \sum^{l - 1}_{k = 0}
  (- 1)^{l + k}  \frac{\q_{k + 1 / 2} - \q_{k - 1 / 2}}{h_k}
\end{aligned}
$$
where $\mathbf{\q}$ is the vector of $\q_l$ for this column and
$\mathbf{d}$ is a vector dependent only on the values of $\q$ in the
neighboring columns. */

double max_slope = 0.577350269189626; // = tan(30.*pi/180.)
#define slope_limited(dz) (fabs(dz) < max_slope ? (dz) :	\
			   ((dz) > 0. ? max_slope : - max_slope))

static void box_matrix (Point point, scalar q, scalar rhs,
			face vector hf, scalar eta,
			double * H, double * d)
{
  coord dz, dzp;
  foreach_dimension()
    dz.x = zb[] - zb[-1], dzp.x = zb[1] - zb[];
  foreach_layer()
    foreach_dimension()
      dz.x += h[] - h[-1], dzp.x += h[1] - h[];
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = h[0,0,m]/(sq(Delta)*cm[]);
    d[l] = rhs[0,0,m];
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    foreach_dimension() {
      double s = Delta*slope_limited((dz.x - h[0,0,m] + h[-1,0,m])/Delta);
      double sp = Delta*slope_limited((dzp.x - h[1,0,m] + h[0,0,m])/Delta);
      d[l] -= a*(gmetric(0)*(h[-1,0,m] - s)*q[-1,0,m] +
		 gmetric(1)*(h[1,0,m] + sp)*q[1,0,m]
#if 1
		 + 2.*G*theta_H*
		 (hf.x[1,0,m]*gmetric(1)*(eta[1,0,m] - eta[0,0,m]) -
		  hf.x[0,0,m]*gmetric(0)*(eta[0,0,m] - eta[-1,0,m]))
#endif
		 );
      H[l*nl + l] -= a*(gmetric(0)*(h[0,0,m] + s) +
			gmetric(1)*(h[0,0,m] - sp));
    }
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      foreach_dimension() {
        double s = Delta*slope_limited(dz.x/Delta);
        double sp = Delta*slope_limited(dzp.x/Delta);
	d[l] -= a*(gmetric(0)*(h[-1,0,m] + s)*q[-1,0,m+1] +
		   gmetric(1)*(h[1,0,m] - sp)*q[1,0,m+1]);
	H[l*(nl + 1) - 1] -= a*(gmetric(0)*(h[0,0,m] - s) +
				gmetric(1)*(h[0,0,m] + sp));
      }
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = h[0,0,nl-1-k];
      if (hk > dry) {
	H[l*nl + k] -= 8.*s*h[0,0,m]/hk;
	H[l*nl + k - 1] += 8.*s*h[0,0,m]/hk;
      }
    }
    foreach_dimension()
      dz.x -= h[0,0,m] - h[-1,0,m], dzp.x -= h[1,0,m] - h[0,0,m];
  }
}

/**
## Relaxation operator

The updated values of $\q$ in a column are obtained as
$$
\mathbf{\q} = \mathbf{H}^{-1}\mathbf{b}
$$
were $\mathbf{H}$ and $\mathbf{b}$ are the Hessenberg matrix and
vector constructed by the function above. */

#include "hessenberg.h"

face vector hf;

#define qflux(qf,i) {							\
  double dz = zb[i] - zb[i-1];						\
  foreach_layer() {							\
    qf = 0.;								\
    if (h[i] + h[i-1] > dry) {						\
      double s = Delta*slope_limited(dz/Delta);				\
      qf += (h[i] + s)*q[i] - (h[i-1] - s)*q[i-1];			\
      if (point.l < nl - 1) {						\
	double s = Delta*slope_limited((dz + h[i] - h[i-1])/Delta);     \
	qf += (h[i] - s)*q[i,0,1] - (h[i-1] + s)*q[i-1,0,1];		\
      }									\
      qf *= gmetric(i)*hf.x[i]/(Delta*(h[i] + h[i-1]));			\
    }

#define end_qflux(i)				\
    dz += h[i] - h[i-1];			\
  }						\
}

trace
static void relax_nh (scalar * ql, scalar * rhsl, int lev, void * data)
{
  scalar q = ql[0], rhs = rhsl[0];
  scalar eta = ql[1], rhs_eta = rhsl[1];
  face vector alpha = *((vector *)data);
  foreach_level_or_leaf (lev) {
    // q
    double H[nl*nl], b[nl];
    box_matrix (point, q, rhs, hf, eta, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
    foreach_layer()
      q[] = b[l--];

    // eta
    double n = 0.;
    foreach_dimension() {
      double qf;
      qflux(qf,0)
	n += qf;
      end_qflux(0);
      qflux(qf,1)
	n -= qf;
      end_qflux(1);
    }
    n *= theta_H*sq(dt)*Delta/cm[];
    n -= sq(Delta)*rhs_eta[];
    double d = - sq(Delta);
    foreach_dimension() {
      n += (alpha.x[1]*eta[1] + alpha.x[]*eta[-1])/cm[];
      d += (alpha.x[1] + alpha.x[])/cm[];
    }
    eta[] = n/d;
  }
}

/**
## Residual computation

The residual is computed as
$$
\begin{aligned}
  \text{res}_l = & \text{rhs}_l -
  h_l \partial_x \partial_x (h_l \q_{l - 1 / 2}) -
  h_l \partial_x \partial_x (h_l \q_{l + 1 / 2}) -\\
  & 4 (\q_{l + 1 / 2} - \q_{l - 1 / 2}) - 8 h_l \sum^{l - 1}_{k = 0}
  (- 1)^{l + k}  \frac{\q_{k + 1 / 2} - \q_{k - 1 / 2}}{h_k}
\end{aligned}
$$
*/

#if 0
trace
static double residual_nh (scalar * ql, scalar * rhsl,
			   scalar * resl, void * data)
{
  scalar q = ql[0], rhs = rhsl[0], res = resl[0];
  scalar eta = ql[1], rhs_eta = rhsl[1], res_eta = resl[1];
  //  face vector alpha = *((vector *)data); // fixme
  double maxres = 0.;
  double C = G*sq(dt);
#if 1 // TREE
  face vector g = new face vector[nl];
  foreach_face() {
    double qf;
    qflux(qf,0) {
      g.x[] = gmetric(0)*(qf + 2.*hf.x[]*theta_H*G*(eta[] - eta[-1]))/Delta;
    } end_qflux(0);
  }
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    /* conservative coarse/fine discretisation (2nd order) */
    // residual for q
    coord dz, dzp;
    foreach_dimension() {
      dz.x = zb[] - zb[-1];
      dzp.x = zb[1] - zb[];
    }
    foreach_layer() {
      res[] = rhs[] + 4.*q[];      
      foreach_dimension() {
	res[] -= h[]*(g.x[1] - g.x[])/(Delta*cm[]);
	res[] += h[]*(g.x[]*slope_limited((dz.x + h[] - h[-1])/Delta) +
		      g.x[1]*slope_limited((dzp.x + h[1] - h[])/Delta))/
	  (hf.x[] + hf.x[1] + dry);
#if 0 // THETA_SCHEME==1 // fixme: why?
	  /theta_H
#endif
	  ;
	if (point.l > 0)
	  res[] -= h[]*(g.x[0,0,-1]*slope_limited(dz.x/Delta) +
			g.x[1,0,-1]*slope_limited(dzp.x/Delta))/
	    (hf.x[0,0,-1] + hf.x[1,0,-1] + dry);
#if 0 // THETA_SCHEME==1 // fixme: why
	    /theta_H
#endif
	    ;
      }
      if (point.l < nl - 1)
        res[] -= 4.*q[0,0,1];
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
	double hk = h[0,0,k];
	if (hk > dry)
	  res[] += 8.*s*(q[0,0,k] - q[0,0,k+1])*h[]/hk;
      }
#if 1
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
#endif
      foreach_dimension() {
	dz.x += h[] - h[-1];
	dzp.x += h[1] - h[];
      }
    }
    // residual for eta
    res_eta[] = rhs_eta[] - eta[];
    foreach_layer()
      foreach_dimension()
        res_eta[] += theta_H*sq(dt)/2.*(g.x[1] - g.x[])/(Delta*cm[]);
#endif // TREE

#if 0    
#if 1 // non-dimensional maximal residual
    if (fabs(res_eta[]/C) > maxres)
      maxres = fabs(res_eta[]/C);
#else
    if (fabs(res_eta[]) > maxres)
      maxres = fabs(res_eta[]);
#endif
#endif
  }
  boundary (resl);
#if 1 // TREE
  delete ((scalar *){g});
#endif
  return maxres;
}
#else
trace
static double residual_nh (scalar * ql, scalar * rhsl,
			   scalar * resl, void * data)
{
  scalar q = ql[0], rhs = rhsl[0], res = resl[0];
  scalar eta = ql[1], rhs_eta = rhsl[1], res_eta = resl[1];
  //  face vector alpha = *((vector *)data); // fixme
  double maxres = 0.;
#if 1 // TREE
  face vector g = new face vector[nl];
  foreach_face() {
    double qf;
    qflux(qf,0) {
      g.x[] = 2.*(qf + gmetric(0)*hf.x[]*theta_H*G*(eta[] - eta[-1])/Delta);
    } end_qflux(0);
  }
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    /* conservative coarse/fine discretisation (2nd order) */
    // residual for q
    coord dz;
    foreach_dimension()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer() {
      res[] = rhs[] + 4.*q[];      
      foreach_dimension() {
	res[] -= h[]*(g.x[1] - g.x[])/(Delta*cm[]);
	res[] += h[]*(g.x[] + g.x[1])/(hf.x[] + hf.x[1] + dry)*
	  slope_limited((dz.x + hf.x[1] - hf.x[])/(Delta*cm[]))
#if 0 // THETA_SCHEME==1 // fixme: why?
	  /theta_H
#endif
	  ;
	if (point.l > 0)
	  res[] -= h[]*(g.x[0,0,-1] + g.x[1,0,-1])/
	    (hf.x[0,0,-1] + hf.x[1,0,-1] + dry)*slope_limited(dz.x/(Delta*cm[]))
#if 0 // THETA_SCHEME==1 // fixme: why
	    /theta_H
#endif
	    ;
      }
      if (point.l < nl - 1)
        res[] -= 4.*q[0,0,1];
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
	double hk = h[0,0,k];
	if (hk > dry)
	  res[] += 8.*s*(q[0,0,k] - q[0,0,k+1])*h[]/hk;
      }
#if 1
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
#endif
      foreach_dimension()
	dz.x += hf.x[1] - hf.x[];
    }
    // residual for eta
    res_eta[] = rhs_eta[] - eta[];
    foreach_layer()
      foreach_dimension()
        res_eta[] += theta_H*sq(dt)/2.*(g.x[1] - g.x[])/(Delta*cm[]);
#endif // TREE

#if 0
#if 1 // non-dimensional maximal residual
    if (fabs(res_eta[]/C) > maxres)
      maxres = fabs(res_eta[]/C);
#else
    if (fabs(res_eta[]) > maxres)
      maxres = fabs(res_eta[]);
#endif
#endif
  }
  boundary (resl);
#if 1 // TREE
  delete ((scalar *){g});
#endif
  return maxres;
}
#endif

event acceleration (i++)
{
  scalar rhs = new scalar[nl];
  double h1 = 0., v1 = 0.;
  foreach (reduction(+:h1) reduction(+:v1)) {
#if 0 // necessary for galilean_invariance.tst
    coord dz;
    foreach_dimension()
      dz.x = zb[1] - zb[-1];
    foreach_layer() {
      rhs[] = 2.*w[];
      foreach_dimension()
	rhs[] -= u.x[]*slope_limited((dz.x + h[1] - h[-1])/(2.*Delta));
      if (point.l > 0)
	foreach_dimension()
	  rhs[] += u.x[0,0,-1]*slope_limited(dz.x/(2.*Delta));
      foreach_dimension()
	rhs[] += ((h[] + h[1])*fm.x[1]*
		  (hf.x[1] > dry ? hu.x[1]/hf.x[1] : 0.) -
		  (h[] + h[-1])*fm.x[]*
		  (hf.x[] > dry ? hu.x[]/hf.x[] : 0.))/(2.*Delta);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
	rhs[] += 4.*s*w[0,0,k];
      rhs[] *= 2.*h[]/(theta_H*dt);
      foreach_dimension()
	dz.x += h[1] - h[-1];
      h1 += dv()*h[];
      v1 += dv();
    }
#elif 1
    coord dz;
    foreach_dimension()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer() {
      rhs[] = 2.*w[];
      // fixme: add comment on Galilean invariance
      // fixme: this should be consistent with the residual above
      foreach_dimension()
	rhs[] -= u.x[]*slope_limited((dz.x + hf.x[1] - hf.x[])/(Delta*cm[]));
      if (point.l > 0)
	foreach_dimension()
	  rhs[] += u.x[0,0,-1]*slope_limited(dz.x/(Delta*cm[]));
      foreach_dimension()
	rhs[] += (hu.x[1] - hu.x[])/(Delta*cm[]);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
	rhs[] += 4.*s*w[0,0,k];
      rhs[] *= 2.*h[]/(theta_H*dt);
      foreach_dimension()
	dz.x += hf.x[1] - hf.x[];
      h1 += dv()*h[];
      v1 += dv();
    }    
#else
    coord dz, dzp;
    foreach_dimension() {
      dz.x = zb[] - zb[-1];
      dzp.x = zb[1] - zb[];
    }
    foreach_layer() {
      rhs[] = 2.*w[];
      foreach_dimension()
	rhs[] -= (hu.x[]*slope_limited((dz.x + h[] - h[-1])/Delta) +
		  hu.x[1]*slope_limited((dzp.x + h[1] - h[])/Delta))/
	(hf.x[] + hf.x[1] + dry);
      if (point.l > 0)
	foreach_dimension()
	  rhs[] += (hu.x[0,0,-1]*slope_limited(dz.x/Delta) +
		    hu.x[1,0,-1]*slope_limited(dzp.x/Delta))/
	  (hf.x[0,0,-1] + hf.x[1,0,-1] + dry);
      foreach_dimension()
	rhs[] += (hu.x[1] - hu.x[])/(Delta*cm[]);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
	rhs[] += 4.*s*w[0,0,k];
      rhs[] *= 2.*h[]/(theta_H*dt);
      foreach_dimension() {
	dz.x += h[] - h[-1];
	dzp.x += h[1] - h[];
      }
      h1 += dv()*h[];
      v1 += dv();
    }
#endif
  }

  scalar res = new scalar[nl];
  mgp = mg_solve ({q,eta}, {rhs,rhs_eta}, residual_nh, relax_nh, &alpha_eta,
		  res = {res,res_eta}, // fixme
		  nrelax = 4, minlevel = 1,
		  tolerance = TOLERANCE*sq(h1/(dt*v1)));
  delete ({rhs, res});

  face vector su[];
  foreach_face() {
    su.x[] = 0.;
    double qf;
    qflux(qf,0) {
      ha.x[] -= qf;
      su.x[] += qf;
      hu.x[] -= theta_H*dt*qf;
    } end_qflux(0);
  }
  boundary_flux ({su});

  foreach() {
    foreach_dimension()
      rhs_eta[] += theta_H*sq(dt)*(su.x[1] - su.x[])/(Delta*cm[]);
    double wmax = 0.;
    foreach_layer()
      wmax += h[];
    wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    foreach_layer()
      if (h[] > dry) {
	if (point.l == nl - 1)
	  w[] += dt*q[]/h[];
	else
	  w[] -= dt*(q[0,0,1] - q[])/h[];
	if (fabs(w[]) > wmax)
	  w[] = (w[] > 0. ? 1. : -1.)*wmax;
      }    
  }
  boundary ({w});
}

/**
## Cleanup

The *w* and *phi* fields are freed. */
      
event cleanup (i = end, last) {
  delete ({w, q});
}

/**
## References

~~~bib
@article{vitousek2013stability,
  title={Stability and consistency of nonhydrostatic free-surface models using the semi-implicit $\theta$-method},
  author={Vitousek, Sean and Fringer, Oliver B},
  journal={International Journal for Numerical Methods in Fluids},
  volume={72},
  number={5},
  pages={550--582},
  year={2013},
  publisher={Wiley Online Library}
}
~~~
*/
