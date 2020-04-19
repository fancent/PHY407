import numpy as np
import matplotlib.pyplot as plt
from random import random

Tmax = 1
Tmin = 1e-3
tau = 1e4
x_0 = 2
y_0 = 2

#def f(x,y):
#	return x**2 - np.cos(4*np.pi*x) + (y-1)**2

def f(x,y):
	if 0<x<50 and -20<y<20:
		return np.cos(x) + np.cos(np.sqrt(2)*x + np.cos(np.sqrt(3)*x)) + (y-1)**2
	else:
		return 1e10
t = 0
T = Tmax
fxy = f(x_0,y_0)
x = x_0
y = y_0

x_results = [x]
y_results = [y]
t_results = [t]

while T>Tmin:
    
    # Cooling
    t += 1
    T = Tmax*np.exp(-t/tau)

    # Choose two cities to swap and make sure they are distinct
    rx = np.random.normal(0,1)
    ry = np.random.normal(0,1)
    oldx = x
    oldy = y
    oldfxy = fxy
    x += rx
    y += ry
    fxy = f(x,y)
    diff = fxy - oldfxy

    t_results.append(t)
    x_results.append(x)
    y_results.append(y)
    
    if random()>np.exp(-diff/T):
        x = oldx
        y = oldy
        fxy = oldfxy

# print("x =", x, "y =", y,"f(x,y) =", fxy)
fig, graph = plt.subplots()
graph.plot(t_results, x_results, ".")
graph.set(xlabel='time', ylabel='x value',
       title="plot of x vs t")
graph.grid()
fig.savefig("q1b_x.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(t_results, y_results, ".")
graph.set(xlabel='time', ylabel='y value',
       title="plot of y vs t")
graph.grid()
fig.savefig("q1b_y.png")
plt.show()
