set size ratio -1
set key out
plot 'cells' w l t '', 'out' w l t "exact", '< grep -v halo log' w l t "VOF", \
     '< grep halo log' u 2:3 w p t 'halo'

