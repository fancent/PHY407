"""
This function aims to simulate Thermal diffusion in the Earth’s crust 
with a given initial conditions. We will be using FTCS method
to calculate the temperatures changes.
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
L = 20        # depth (m)
D = 0.1       # m^2 day^-1
A = 10        # Temperature constant (C)
B = 12        # Temperature coeffiecient (C)
N = 100       # Number of divisions in grid
a = L/N       # Grid spacing
h = 0.01      # Time-step
year = 365    # days in a year
c = h*D/(a*a) # FTCS coeffiecient
depth=np.arange(0,20+a,a)

def meanTemp(t):
    return A+(B*np.sin(2*np.pi*t/year))

# Create the arrays
T = np.empty(N+1,float)
T[N] = 11           # constant temperature at 20m
T[0] = meanTemp(0)   # temperature at the surface at t = 0 day
T[1:N] = 10         # initial temperature everywhere else

Tp = np.empty(N+1,float)
Tp[N] = 11           # constant temperature at 20m
Tp[0] = meanTemp(0)   # temperature at the surface at t = 0 day

# Create time intervals for part b
march = 9*year + year//4
june = 9*year + 2 * year//4
sep = 9*year + 3 * year//4
dec = 10*year

def FTCS(T, Tp, start, end):
    """
    This function basically performs FTCS method using the T and Tp arrays
    and updating from start time to end time.
    """
    time = start
    while time < end:
        # Update or maintain end points values 
        T[0] = meanTemp(time)
        T[N] = 11
        Tp[1:N] = T[1:N] + c*(T[0:N-1] + T[2:N+1] - 2 * T[1:N])
        T,Tp = Tp,T
        time += h
    return T, Tp

# PART A
nineYears, Tp1 = FTCS(T,Tp, 0, 9*year)
# plotting
# NOTE: I adjusted the xtick to match with how the array is set up
fig, graph = plt.subplots()
graph.plot(depth,nineYears)
graph.set(xlabel='Depth (m)', ylabel='Temperature (C)',
       title="Temperature at difference depth\nof Earth\'s crust after 9 years")#,
#       xticklabels=['', 0, 4, 8, 12, 16, 20])
graph.grid()
fig.savefig("q2a.png")
plt.show()

# PART B
fig, graph = plt.subplots()

# calculating temperature profile every 3 months and plotting
tenth_march, Tp2 = FTCS(nineYears, Tp1, 9*year, march)
graph.plot(depth,tenth_march, label="1st quarter")

tenth_june, Tp3 = FTCS(tenth_march, Tp2, march, june)
graph.plot(depth,tenth_june, label="2nd quarter")

tenth_sep, Tp4 = FTCS(tenth_june, Tp3, june, sep)
graph.plot(depth,tenth_sep, label="3rd quarter")

tenth_dec, Tp5 = FTCS(tenth_sep, Tp4, sep, dec)
graph.plot(depth,tenth_dec, label="4th quarter")

graph.set(xlabel='Depth (m)', ylabel='Temperature (C)',
       title="Temperature at difference depth of Earth\'s\ncrust in quarters during the 10th year"
       )
graph.grid()
graph.legend()
fig.savefig("q2b.png")
plt.show()