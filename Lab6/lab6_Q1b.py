"""
We want to use a similar code to part a, but modify the function
used at the beginning to give a van der Pol oscillation simulation through our plots
this time. We also want to vary the mu value in this new function
for 3 different values, and then plot the results along side our original 
functions for reference
"""
import numpy as np
import matplotlib.pyplot as plt
#Question 1b
w = 1.0
mu = 0.1 #define the new constant
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    #Modified van der Pol dvelocity equation
    fv = mu*(1-x**2)*v-(w**2)*x
    return np.array([fx,fv],float)
#Including the code from part a) to include in plot for reference
def g(r1,t):
    x = r1[0]
    v = r1[1]
    gx = v
    gv = -(w**2)*x
    return np.array([gx,gv],float)
        
    
a = 0.0
b = 60.0
N = 1000
h = (b-a)/N

#Create distinct position and velocity bins for regular SOH and van der Pol
time = np.arange(a,b,h)
#Van der Pol v and x
xvals = [] 
vvals = []
#SOH v and x
x1vals = []
v1vals = []

#Van der Pol RK4 ODE  solution as in part a
r = np.array([1.0,0.0],float)
for t in time:
    xvals.append(r[0])
    vvals.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+ 0.5*h)
    k3 = h*f(r+0.5*k2, t +0.5*h)
    k4 = h*f(r + k3, t + h)
    r += (k1+2*k2+2*k3+k4)/6
   
#SOH RK$ ODE solution as in part a
r1 = np.array([1.0,0.0],float)
for t in time:
    x1vals.append(r1[0])
    v1vals.append(r1[1])
    k1 = h*g(r1,t)
    k2 = h*g(r1+0.5*k1,t+ 0.5*h)
    k3 = h*g(r1+0.5*k2, t +0.5*h)
    k4 = h*g(r1 + k3, t + h)
    r1 += (k1+2*k2+2*k3+k4)/6   
 
#Plotting both functions
fig, graph = plt.subplots()
graph.plot(time, xvals, label = "Van der Pol")
graph.plot(time,x1vals, label = "SOH")
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'position vs time, van der Pol mu = 0.1')
graph.legend()
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals, label = "Van der Pol")
plt.plot(x1vals, v1vals, label = "SOH")
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'phase space plot, van der Pol mu = 0.1')
graph.legend()
plt.show()

"""
from here on, the code is repeated for the next two desired values of mu, the only change is the constant
mu defined at the beginning, the rest is exactly as above
"""
#mu = 2
w = 1.0
mu = 2.0
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    fv = mu*(1-x**2)*v-(w**2)*x
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
graph.plot(time, xvals, label = "Van der Pol")
graph.plot(time,x1vals, label = "SOH")
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'position vs time, van der Pol mu = 2')
graph.legend()
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals, label = "Van der Pol")
plt.plot(x1vals, v1vals, label = "SOH")
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'phase space plot, van der Pol mu = 2')
graph.legend()
plt.show()

#mu = 5
w = 1.0
mu = 5.0
def f(r,t):
    x = r[0]
    v = r[1]
    fx = v
    fv = mu*(1-x**2)*v-(w**2)*x
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
graph.plot(time, xvals, label = "Van der Pol")
graph.plot(time,x1vals, label = "SOH")
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'position vs time, van der Pol mu = 5')
graph.legend()
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals, label = "Van der Pol")
plt.plot(x1vals, v1vals, label = "SOH")
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'phase space plot, van der Pol mu = 5')
graph.legend()
plt.show()
