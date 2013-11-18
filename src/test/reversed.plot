reset
set title 'Time-reversed advection'

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
plot 'log' u 1:4 t 'max', 'log' u 1:2 t 'norm1', exp(f(log(x))) t ftitle(a,b), exp(f2(log(x))) t ftitle(a2,b2)

if (batch) set term @PNG; set output "interface.png"; else pause -1;
reset
set title 'Time-reversed advection'
set size ratio -1
plot [-0.5:0.5][-0.5:0.5]'out' w l t ""
