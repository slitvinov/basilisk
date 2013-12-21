Amp = 0.01
e2 = 0.1
set yrange [0.008:]
set ylabel 'Q'
set xlabel 'x'
plot 'log' w l t 'numerical', Amp*exp(-e2/2.*x) t 'linear theory';
