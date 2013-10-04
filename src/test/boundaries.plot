set size ratio -1
set key outside top right
set yrange [-0.1:1.1]
plot 'out' w l t '', \
     '< grep bottom log' u 2:3:(5-$4) ps variable pt 4 t 'bottom', \
     '< grep right log' u 2:3:(5-$4) ps variable pt 6 t 'right', \
     '< grep post log' u 2:3:(5-$4) ps variable pt 12 t 'post', \
     '< grep halo log' u 2:3:(5-$4) ps variable pt 14 t 'halo'
