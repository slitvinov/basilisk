set xlabel 'y'
plot 'log' u 2:3 t "", \
     -0.51*(x - 0.5) + x**2/2. - 0.125 t "hydrostatic pressure", \
     0.01 + 0.5 - x t "density"
