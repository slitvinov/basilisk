set xrange [0:1]
set xlabel 'r'

set ylabel 'rho'
plot './out' u 1:2 w p pt 7 ps 0.2 t 'Adaptive', \
     './cout' u 1:2 w p pt 7 ps 0.2 t 'Cartesian'

set output 'velocity.png'
set ylabel 'Normal velocity'
plot './out' u 1:3 w p pt 7 ps 0.2 t 'Adaptive', \
     './cout' u 1:3 w p pt 7 ps 0.2 t 'Cartesian'
