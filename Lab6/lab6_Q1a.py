"""
We want to create a plot of a simple harmonic oscilator showing position vs time, and
velocity vs position with 
a set frequency over a time of 60 seconds starting from a stationary point at x=1
using the RK4 method of differential equation estimation
"""
import numpy as np
import matplotlib.pyplot as plt
#Question 1a
w = 1.0 #omega
#Create a function for our SOH
def f(r,t): 
    #we create two columns within "r" so that we include both position and velocity in our calculations
    x = r[0]
    v = r[1]
    fx = v
    fv = -(w**2)*x
    return np.array([fx,fv],float)
    
#Define constants
a = 0.0
b = 60.0
N = 1000
h = (b-a)/N

#create bins for our position, time and velocity values
time = np.arange(a,b,h)
xvals = []
vvals = []

#initial conditions for our ODE estimate using the RK4 approximation
r = np.array([1.0,0.0],float)
# RK4 code using form derived from text; Mark Newman Computational Physics page 350
#Appending the solved values to the x and v bins
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
graph.set(xlabel = 'time (s)', ylabel = 'position (m)', title = 'omega = 1 rad/s, position vs time, SOH')
plt.show()
fig, graph = plt.subplots()
plt.plot(xvals, vvals)
graph.set(ylabel = 'velocity (m/s)', xlabel = 'position (m)', title = 'omega = 1 rad/s, phase space plot, SOH')
plt.show()
