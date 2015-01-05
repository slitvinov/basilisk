/**
# Using discharge function

In this example, we use the discharge() function in order to impose
different inflows on different rivers situated on the same border. We
include the Saint venant solver and the file discharge.h . We also
define River1[] and River2[] to define the two differents river beds.  */

/**  
# Declarations ... 
*/
#include "saint-venant.h"
#include "discharge.h"

#define LEVEL 7

double pause,a1,a2; 
scalar River1[],River2[];
int imax;

/**
# Parameters

Definition of parameters and calling of the saint venant subroutine run().
*/
int main()
{  
  L0 = 10.;
  X0 = -L0/2.;
  Y0 = -L0/2.;
  G = 9.81;
  N = 1 << LEVEL; 
  pause = 0.0;
  imax=800;
  run();
}

/**
# Initial conditions

We define a topography with two different river beds.  We also define the
scalars River1 and River2 which are equal to 1 on their own bed and to
zero otherwise. We must add a little volume of water near the boundary
in order to impose a time scale to the solver */
event init (i = 0)
{
  double dx=L0/((double)N);
  foreach(){
    zb[] = 0.05*x*x*x*x-x*x+2+0.2*(y+Y0);
    u.x[]=0;
    u.y[]=0;
    h[]= ( y>=Y0+L0-dx ) ? max(-1-zb[],0) : 0;

    // The river n°1 is on the left of our domain (x<0)
    River1[]= x<0 ? 1:0;
    // And the river 2 on the right (x>0)
    River2[]= x>0 ? 1:0;
  }
  boundary(all);
}


/**
# Boundary conditions :

We impose no condition on the bottom boundary (free outflow
condition). At the top boundary, we impose a Neumann's condition equal
to zero on the y component of the velocity and we fix the x component
of the velocity to zero. The condition on the water height will be
fixed just after by the discharge() function.  */

h[bottom]=h[];
u.x[bottom]=u.x[];
u.y[bottom]=u.y[];

u.y[top]= neumann(0);
u.x[top]= dirichlet(0);

/**
# Inflow

This event calculates the imposed height on each river to have the
good inflow. The inflow of the river n°1 starts to zero and linearly
increases to reach a value of 4 cumsec (m³/s) at 0.2 seconds. The
inflow of the river n°2 starts to zero and linearly increases to reach
a value of 2 cumsec at 0.2 seconds */

event Discharge (i++){
  double Q1 = min ( 4 * t / 0.2, 4 );
  double Q2 = min ( 2 * t / 0.2, 2 );

  a1 = discharge ( Q1, top, River1 );
  a2 = discharge ( Q2, top, River2 );

  h[top] = dirichlet ( max ( a1*River1[] + a2*River2[] - zb[], 0 ) );
}

/**
# Volume output  

We write the total volume of fluids, the volume of the first river and
the one in the second in the error file (log) */

event Volume (i++){
  double volumetot = 0, volume1 = 0 , volume2 = 0;
  double dx = L0/((double)N);  
  foreach(){
    if(x<0) volume1 += h[]*dx*dx;
    else volume2 += h[]*dx*dx;
    volumetot += h[]*dx*dx;
  }
  fprintf(stderr,"%g  %g  %g  %g \n",t,volumetot,volume1,volume2);
}

/**
# Gnuplot output
*/
event initplot(i = 0) {
  printf("set view 80,05\n"
	 "set xlabel 'X'\n"
	 "set ylabel 'Y'\n"
	 "set zlabel 'Hauteur'\n"
	 "set hidden3d; unset ytics ; unset xtics\n");
}
event plot (i<=imax;i+=10 ) {
  double dx=2.*L0/pow(2.,LEVEL),dy=dx;
  printf("set title ' --- t= %.3g '\n"
	 "sp[%g:%g][%g:%g][-5:5]  '-' u 1:2:($3+$4-.05) t'free surface' w l lt 3,'' u 1:2:4 t'topo' w l lt 2\n",t,X0,-X0,Y0,-Y0);
  for(double x=X0;  x<=X0+L0;x+=dx)
    {for(double y=Y0;y<=Y0+L0;y+=dy)
	{ 
	  printf (" %g %g %g  %g \n", x, y, interpolate (h, x, y),  interpolate (zb, x, y) );}
      printf ("\n");
    } 
  printf("e\n"
	 "pause %.5lf \n\n",pause);
}


/**
# Results 
~~~gnuplot Movie
set term gif animate
set output 'Movie.gif'
load './out';

~~~

~~~gnuplot Total volume of fluid
set title 'Total'
set xlabel 'T' 
set ylabel 'Volume' 
set xtics; set ytics
set terminal png
set output 'total.png'
 
f(x)=a*x+b 
fit [0.2:] f(x) './log' u 1:2 via a,b 
tit=sprintf("best fit with a = %1.3f",a) 
plot './log' u 1:2 w p pt 7 ps 0.2 t 'tot vol', \
     f(x) t tit;

~~~

~~~gnuplot River n°1
set title 'River 1'
set xlabel 'T'
set ylabel 'Volume'
set xtics; set ytics
set terminal png
set output 'river1.png'

f(x)=a*x+b 
fit [0.2:] f(x) './log' u 1:3 via a,b
tit=sprintf("best fit with a = %1.3f",a)
plot './log' u 1:3 w p pt 7 ps 0.2 t 'river1 vol', \
    f(x) t tit;

~~~



~~~gnuplot River n°2
set title 'River 2'
set xlabel 'T'
set ylabel 'Volume'
set xtics; set ytics;
set terminal png
set output 'river2.png'

f(x)=a*x+b 
fit [0.2:] f(x) './log' u 1:4 via a,b
tit=sprintf("best fit with a = %1.3f",a)
plot './log' u 1:4 w p pt 7 ps 0.2 t 'river2 vol', \
     f(x) t tit;

~~~

*/
