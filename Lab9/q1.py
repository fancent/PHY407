"""
The purpose of this code is to simulate shallow
water system using the Lax-Wendroff scheme.
Comment and uncomment for both part b and c
"""
import numpy as np
import matplotlib.pyplot as plt

#constants
g = 9.81 
xmin = 0.0
xmax = 1.0
t = 0.0
tmax = 4.0

# for part (b)
dx = 0.02
dt = 0.01

# for part (c)
# dx = 1./150
# dt = 0.001

# Initialize the arrays
col_len = int(xmax/dx)+1
u = np.zeros((col_len, 2))
F = np.zeros((col_len, 2))
xgrid = np.arange(0,xmax+dx, dx)

# for part (b)
H = [0.01 for i in range(len(xgrid))]
# for part (c)
# H = 0.001 + 0.1*np.exp(-7*xgrid)

#Boundary and initial condition
# for part (b)
A = 0.002
mu = 0.5
sigma = 0.05
# for part (c)
# A = 0.0005
# mu = 0.0
# sigma = 0.1
u[:,1] = A*np.exp(-((xgrid-mu)**2)/(sigma**2))

def F_vector(u, eta, h):
    """
    This function calculates the F vector with
    u and eta using H defined above.
    """
    first = 0.5*u**2 + g*eta
    second = (h+ eta)*u
    return first, second

while t <= tmax:
    """
    This while loop updates both u and F vector from initial t to desired end time.
    """
    #forward difference
    F[:,0], F[:,1] = F_vector(u[:,0], u[:,1], H)
    
    #calculating all halfsteps for u vector
    u_j_plus_half = np.zeros((col_len-1,2))
    for j in range(len(u_j_plus_half)):
        u_j_plus_half[j] = 0.5*(u[j+1] + u[j]) - (dt/(2*dx))*(F[j+1]-F[j])

    #calculating all halfsteps for F vector
    F_half = np.zeros((col_len-1,2))
    for i in range(len(u_j_plus_half)):
        F_half[i,0], F_half[i,1] = F_vector(u_j_plus_half[i,0], u_j_plus_half[i,1], H[i])

    #finding n+1 values with the halfsteps
    for k in range(1, len(u)-1):
        u[k] = u[k] - (dt/dx)*(F_half[k]-F_half[k-1])
    #boundary calculations
    u[0,1] = u[0,1] - (dt/dx)*(F[1,1]-F[0,1])
    u[-1,1] = u[-1,1]- (dt/dx)*(F[-1,1]-F[-2,1])

    # for different desired time for part b and c
    # if np.abs(t) < 1*dt or np.abs(t-1.0)< 1.*dt or np.abs(t-4.0) < 1*dt:
    if np.abs(t) < 1*dt or np.abs(t-1.0)< 1.*dt or np.abs(t-2.0)< 1.*dt or np.abs(t-3.0)< 1.*dt or np.abs(t-4.0) < 1*dt:
        plt.plot(xgrid, u[:,1])
        plt.title("η at t={}s".format(t))
        plt.xlabel("Domain (m)")
        plt.ylabel("η / height (m)")
        plt.show()

    t = round(t+dt,3)

