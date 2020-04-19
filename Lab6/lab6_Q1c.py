"""
Here we want to explore the changes in our outputs by incorporating a sinusoidal driver in our oscillator.
The only changes we want to make here compared to the last two parts are;
the addition of a sinusoidal driving term in our initial function for the van der Pol oscillator
and a new constant A (amplitude). Additionally we want to incorporate another omega term
for the sine function in the driver. We will hold mu constant at 1.0, and vary the driving omega term
for three sets of conditions.
"""
import numpy as np
import matplotlib.pyplot as plt
from math import sin
#Question 1c
w = 1.0
w2 = 4.0*w #defining new omega term with first desired condition
mu = 1.0
A = 8.53 #New amplitude constant
#We adjust the dvelocity function for our new driving term
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    fv = A*sin(w2*t)+mu*(1-x**2)*v-(w**2)*x
    return np.array([fx,fv],float)
    
# everything from here on is taken from part a
a = 0.0
b = 60.0
N = 1000
h = (b-a)/N

time = np.arange(a,b,h)
xvals = []
vvals = []

r = np.array([1.0,0.0],float)
for t in time:
    xvals.append(r[0])
    vvals.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+ 0.5*h)
    k3 = h*f(r+0.5*k2, t +0.5*h)
    k4 = h*f(r + k3, t + h)
    r += (k1+2*k2+2*k3+k4)/6
   
 #plotting
fig, graph = plt.subplots()
graph.plot(time, xvals)
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'w2 = 4.0w, position vs time, van der Pol with sinusoidal driver')
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals)
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'w2 = 4.0w, phase space plot, van der Pol with sinusoidal driver')
plt.show()

"""
We repeat the exact code as that above two more times, with two different values of omega2
"""
w = 1.0
w2 = 0.5*w
mu = 1.0
A = 8.53
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    fv = A*sin(w2*t)+mu*(1-x**2)*v-(w**2)*x
    return np.array([fx,fv],float)
    
    
a = 0.0
b = 60.0
N = 1000
h = (b-a)/N

time = np.arange(a,b,h)
xvals = []
vvals = []

r = np.array([1.0,0.0],float)
for t in time:
    xvals.append(r[0])
    vvals.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+ 0.5*h)
    k3 = h*f(r+0.5*k2, t +0.5*h)
    k4 = h*f(r + k3, t + h)
    r += (k1+2*k2+2*k3+k4)/6
   
fig, graph = plt.subplots()
graph.plot(time, xvals)
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'w2 = 0.5w, position vs time, van der Pol with sinusoidal driver')
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals)
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'w2 = 0.5w, phase space plot, van der Pol with sinusoidal driver')
plt.show()

w = 1.0
w2 = 1.0*w
mu = 1.0
A = 8.53
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    fv = A*sin(w2*t)+mu*(1-x**2)*v-(w**2)*x
    return np.array([fx,fv],float)
    
    
a = 0.0
b = 60.0
N = 1000
h = (b-a)/N

time = np.arange(a,b,h)
xvals = []
vvals = []

r = np.array([1.0,0.0],float)
for t in time:
    xvals.append(r[0])
    vvals.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+ 0.5*h)
    k3 = h*f(r+0.5*k2, t +0.5*h)
    k4 = h*f(r + k3, t + h)
    r += (k1+2*k2+2*k3+k4)/6
   
fig, graph = plt.subplots()
graph.plot(time, xvals)
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'w2 = w, position vs time, van der Pol with sinusoidal driver')
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals)
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'w2 = w, phase space plot, van der Pol with sinusoidal driver')
plt.show()