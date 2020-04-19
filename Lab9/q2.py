"""
The Purpose of this program is to use 3 different methods in order to solve FODEs associated with 
a piano string. The three methods are the FTCS, the Crank-Nicolson method and the spectral method. 
We also waant to plot each of these solutions for 5 distinct time values and compare the effectiveness
of each method
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, irfft

#Constants
h = 1e-6 #time step
L = 1 #length of string
v = 100 #velocity m/s
d = 0.1 #m
C = 1 #m/s
s = 0.3 #m
N = 100 
a = L/N
epsilon = h/1000 #error

#defining our velocity function
def psi(x):
	return C*x*(L-x)/L**2*np.exp(-(x-d)**2/2/s**2)
"""
Question 2b
"""


def FTCS(phi,dphi):
    """
    This function applies the FTCS method of solving FODEs for a given
    wave function with specified initial conditions for position and velocity
    """
    phi_old = phi.copy()
    phi[1:N] = phi[1:N] + h*dphi[1:N]
    dphi[1:N] = dphi[1:N] + h*v**2/(a**2)*(phi_old[2:N+1] + phi_old[0:N-1] - 2*phi_old[1:N])
    return phi , dphi

def plotting(phi, x, t):
    """
    This function is meant to help plot our sting waves at the desired time instances
    while helping with a rounding issue in python
    """
    plt.plot(x, phi)
    plt.title("graph at t={}s".format(round(t,4)))
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    plt.show()

phi_0 = np.zeros(N+1,float) #initial conditions for phi
x = np.arange(0,L+a,a)
dphi_0 = psi(x) #derivative of phi is represented by psi

def FTCS_method(t, tmax, phi, dphi):
    """
    This function is meant to apply the FTCS function created earlier, and plot
    the desired time values on the same graph, as well as seperate graphs for comparison
    """
    phi_plots = []
    while t <= tmax:
        phi,dphi = FTCS(phi, dphi)
        
      #Desired time values
        t1 = 0.002
        t2 = 0.004
        t3 = 0.006
        t4 = 0.012
        t5 = 0.1
        
      #Python doesnt round these numbers off, so we allow an error when selecting our values
        if abs(t - t1)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t2)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t3)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t4)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        # if abs(t - t5)<epsilon:
            # plotting(phi, x, t)

        t+=h
    for i in range(len(phi_plots)):
        plt.plot(x, phi_plots[i][1], label="at t={}s".format(round(phi_plots[i][0],4)))
    plt.title("Phi at all t using FTCS method")
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    plt.legend()
    plt.show()

# FTCS_method(0, 0.11, phi_0, dphi_0)

"""
Question 2 c
"""
# creating diagonal matrices for Crank–Nicolson method
Nrange = np.arange(N-2)
coefficientMatrix = -2 * np.identity(N-1)
coefficientMatrix[Nrange, Nrange+1] = 1
coefficientMatrix[Nrange+1, Nrange] = 1

# results derived from substitute equations and rearranging
A = (2/h)*np.identity(N-1) - (h*v**2/(2*a**2)*coefficientMatrix)
A_inv = np.linalg.inv(A)
B = (2/h)*np.identity(N-1) + (h*v**2/(2*a**2)*coefficientMatrix)


def Crank(phi, psi):
    """
    This function calculates the new phi with Crank–Nicolson method
    and then using the new phi to update psi at next time step.
    """
    phi_old = phi.copy()
    phi[1:-1] = np.dot(np.dot(A_inv,B),phi[1:-1]) + 2*np.dot(A_inv, psi[1:-1])
    psi = (2/h)*(phi-phi_old)-psi
    return phi, psi

def Crank_method(t, tmax, phi, dphi):
    """
    This function performs Crank–Nicolson method over a period of time.
    """
    phi_plots = []
    while t <= tmax:
        phi,dphi = Crank(phi, dphi)
        
        t1 = 0.002
        t2 = 0.004
        t3 = 0.006
        t4 = 0.012
        t5 = 0.1
        
        if abs(t - t1)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t2)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t3)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t4)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
        if abs(t - t5)<epsilon:
            phi_plots.append((t, phi.copy()))
            # plotting(phi, x, t)
    
        t+=h
    for i in range(len(phi_plots)):
        plt.plot(x, phi_plots[i][1], label="at t={}s".format(round(phi_plots[i][0],4)))
    plt.title("Phi at all t using Crank–Nicolson method")
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    plt.legend()
    plt.show()
    
# Crank_method(0, 0.11, phi_0, dphi_0)


"""
Question 2 e
"""
#copy the dst and idst functions from dcst.py
def dst(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[0] = y2[N] = 0.0
    y2[1:N] = y[1:]
    y2[:N:-1] = -y[1:]
    a = -np.imag(rfft(y2))[:N]
    a[0] = 0.0

    return a

def idst(a):
    N = len(a)
    c = np.empty(N+1,complex)
    c[0] = c[N] = 0.0
    c[1:N] = -1j*a[1:]
    y = irfft(c)[:N]
    y[0] = 0.0

    return y

#Define k necessary for spectral method, count to N+1 (include bounds)
k = np.arange(N+1)

#Define angular frequency using value from text
omega = np.pi*v*k/L

#first we want to discrete sine transform our velocity to use it to calculate
# our coefficients
dstpsi = dst(dphi_0)

#Get solutions by taking the inverse of the coefficient
def spectrum(t):
    #We ignore the first chunk of the coefficient because it goes to zero as the initial phi is zero
    phi = idst((dstpsi/omega)*np.sin(omega*t)) 
    return phi
    
#Here we plot our solutions the same way as we did previously
def spectrum_method(t,tmax):
    phi_plots = []
    while t <= tmax:
        phi = spectrum(t)
        
        t1 = 0.002
        t2 = 0.004
        t3 = 0.006
        t4 = 0.012
        t5 = 0.1
        
        if abs(t - t1)<epsilon:
            phi_plots.append((t, phi))
        if abs(t - t2)<epsilon:
            phi_plots.append((t, phi))   
        if abs(t - t3)<epsilon:
            phi_plots.append((t, phi))  
        if abs(t - t4)<epsilon:
            phi_plots.append((t, phi))   
        if abs(t - t5)<epsilon:
            phi_plots.append((t, phi))
        

        t+=h
    for i in range(len(phi_plots)):
        plt.plot(x, phi_plots[i][1], label="at t={}s".format(round(phi_plots[i][0],5)))
    plt.title("Phi at all t using spectral method")
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    plt.legend()
    plt.show()
    
spectrum_method(0, 0.11)