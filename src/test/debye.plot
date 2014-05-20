set xlabel 'x'
gamma = tanh(0.25)
fi(x) = 2*log((1+gamma*exp(-sqrt(2)*x))/(1-gamma*exp(-sqrt(2)*x)))
nplus(x) = exp(-fi(x))
nminus(x) = exp(fi(x))
plot 'log' u 1:2 notitle, fi(x) t '{/Symbol f}',\
     'log' u 1:3 notitle, nplus(x) t 'n+',\
     'log' u 1:4 notitle, nminus(x) t 'n-'