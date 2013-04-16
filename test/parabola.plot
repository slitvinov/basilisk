reset
set title 'Oscillations in a parabolic container'

h0 = 10.
a = 3000.
tau = 1e-3
B = 5.
G = 9.81
p = sqrt (8.*G*h0)/a
s = sqrt (p*p - tau*tau)/2.
u0(t) = B*exp (-tau*t/2.)*sin (s*t)

set xlabel 'x (m)'
set ylabel 'z (m)'
t = 1500
psi(x) = a*a*B*B*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + \
      (tau*tau/4. - s*s)*cos (2.*s*t)) - B*B*exp(-tau*t)/(4.*G) - \
      exp (-tau*t/2.)/G*(B*s*cos (s*t) + tau*B/2.*sin (s*t))*x + h0
bed(x) = h0*(x/a)**2
set key top center
plot [-5000:5000] \
      '< grep ^p parabola.out' u 2:5:($5+$3) w filledcu lc 3 t 'Numerical', \
      psi(x) > bed(x) ? psi(x) : bed(x) lc 2 t 'Analytical', \
      bed(x) lw 3 lc 1 lt 1 t 'Bed profile'

if (batch) set term pngcairo; set output "parabola_u0.png"; else pause -1;
reset
set key top right
set ylabel 'u0'
set xlabel 'Time'
plot u0(x) t 'Analytical', '< grep ^s parabola.out' u 2:(-$3) every 2 w p t 'Numerical'

if (batch) set term pngcairo enhanced; set output "parabola_convergence.png"; else pause -1;
reset
set xlabel 'Resolution'
set ylabel 'Relative error norms'
set key bottom left
set logscale
set xtics 32,2,512
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f1(x)=a1+b1*x
fit f1(x) 'parabola.log' u (log($1)):(log($2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) 'parabola.log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'parabola.log' u (log($1)):(log($4)) via am,bm
plot exp (f1(log(x))) t ftitle(a1,b1), \
     exp (f2(log(x))) t ftitle(a2,b2), \
     exp (fm(log(x))) t ftitle(am,bm),  \
     'parabola.log' u 1:2 t '|h|_1' ps 1.5, \
     'parabola.log' u 1:3 t '|h|_2' ps 1.5, \
     'parabola.log' u 1:4 t '|h|_{max}' ps 1.5 lc 0
