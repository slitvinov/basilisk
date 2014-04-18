import matplotlib.pyplot as plt
import numpy as np
import stream as bas

def init(i,t):
    bas.omega.f = bas.noise

N = 256
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)

def graph(i,t):
    print "t=",t
    Z = bas.omega.f(X,Y)
    plt.imshow(Z)
    plt.draw()
    plt.show(block=False)

bas.resolution(N)
bas.event(init, i = 0)
bas.event(graph, t = range(0,1000,10))
bas.run()
plt.show()
