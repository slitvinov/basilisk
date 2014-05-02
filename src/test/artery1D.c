/**
# 1D arterial flow

A 1D model for arterial flows can be derived from the Navier-Stokes
equations, in terms of the cross sectional area $A$ and flow rate $Q$,
we have
$$
\partial_t A +\partial_x Q  = 0 
$$
$$
\partial_t Q +\partial_x (Q^2/A) = - A \partial_x p/\rho - f_r 
$$ 
where $p(A)$ models the wall properties of the arteries, $\rho$ is the
blood density and $f_r$ stands for the wall shear stress. For a simple
linear wall relation, $p = K A$ with $K$ a constant, we can write the
flux as $F = (Q,Q^2/A + 2 e_1 A)$ and the source term as $S = (0,-e_2
Q/A)$ using two parameters $e_1$ and $e_2$.

We will use a 1D Cartesian grid and the generic solver for systems of
conservation laws. */

#include "grid/cartesian1D.h"
#include "conservation.h"

/**
## Variables

We define the conserved scalar fields $a$ and $q$ which are passed to the
generic solver through the *scalars* list. We don't have any conserved
vector field. */

scalar a[], q[];
scalar * scalars = {a,q};
vector * vectors = NULL;

/**
The other parameters are specific to the example. */

double e1, e2, omega, Amp;

/**
## Functions

We define the *flux* function required by the [generic
solver](/src/conservation.h). */

void flux (const double * s, double * f, double e[2])
{  
  double a = s[0], q = s[1], u = q/a;
  f[0] = q;
  f[1] = q*q/a + e1*a*a;
  // min/max eigenvalues
  double c = sqrt(2.*e1*a);
  e[0] = u - c; // min
  e[1] = u + c; // max
}

/**
We need to add the source term of the momentum equation. We first
define a function which, given the current states, fills the *updates*
with the source terms for each conserved quantity. */

static void momentum_source (scalar * current, scalar * updates)
{

  /**
  We recover the current fields and their variations from the lists... */

  scalar a = current[0], q = current[1], da = updates[0], dq = updates[1];

  /**
  We initialise the source terms for each conserved variable. Note that
  *a* and *da*, *q* and *dq* may be the same fields, so that care needs
  to be taken with the order of operations. */

  foreach() {
    dq[] = - e2*q[]/a[];
    da[] = 0.;
  }
}

/**
We then overload the default *sources* function pointer of the generic
solver. */

event defaults (i = 0)
  sources = momentum_source;

/**
## Boundary conditions

We impose a sinusoidal flux $Q(t)$ at the left of the domain. */

q[left] = dirichlet(Amp*sin(2.*pi*omega*t));

/**
## Parameters

For small amplitudes $Amp = 0.01$ at the input boundary condition the
system has analytical solutions for $e1 < e2$, in this case the spatial
envelope of the flux rate behaves like $Q=Amp\times e^{-e2/2x}$ [Wang at
al., 2013]. */

int main() {
  init_grid (400);
  e1 = 0.5 ;
  e2 = 0.1 ;
  omega = 1.;
  Amp = 0.01 ;
  run();
}

/**
## Initial conditions 

The initial conditions are $A=1$ and $Q=0$. */

event init (i = 0) {
  foreach()
    a[] = 1.;
}

/**
## Outputs

We print to standard error the spatial profile of the flow rate
$Q$. */

event printdata (t += 0.1; t <= 1.) {
  foreach()
    fprintf (stderr, "%g %.6f \n", x, q[]);
  fprintf (stderr, "\n\n");
}

/**
We get the following comparison between the numerical solution and the
linear theory for the flow rate $Q$.

~~~gnuplot
set output 'plot.png'
Amp = 0.01
e2 = 0.1
set yrange [0.008:]
set ylabel 'Q'
set xlabel 'x'
plot 'log' w l t 'numerical', Amp*exp(-e2/2.*x) t 'linear theory'
~~~
*/

