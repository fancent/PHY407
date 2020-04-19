import numpy as np
from scipy import constants
from scipy import special
import matplotlib.pyplot as plt

Q = 10**-13
l = 1

#defining simpsons rule
def simpsonsRule(f,a,b,N):
    h = (b-a)/N
    oddSum = 0
    for k in range(1, N, 2):
        oddSum += f(a+k*h)
    evenSum = 0
    for k in range(2, N, 2):
        evenSum += f(a+k*h)
    integral = (h/3)*(f(a)+f(b)+4*oddSum+2*evenSum)
    return integral

#defining constants for simpson's rule where N is the number of slices
#a is the lower bound and b is the upper bound
N = 60
a = -constants.pi/2
b = constants.pi/2

#defining eq 3
def V(z, r):
    #function generator for each z and r
    def f(u):
        return (Q*np.exp(-(np.tan(u))**2))/(4*np.pi*constants.epsilon_0*(np.cos(u)**2)*np.sqrt(((z-l*np.tan(u))**2)+r**2))
    #calculate the integral
    integral = simpsonsRule(f, a, b, N)
    return integral

#defining eq 4
def Vzero(r):
    return (Q/(4*np.pi*constants.epsilon_0*l))*np.exp(r**2/(2*l**2))*special.kn(0,r**2/(2*l**2))

#initialize lists to store respective results
VResults = []
VZeroResults = []
#creating a data set of radius from 0.25mm to 5.0mm with 238 entries
radiusRange = np.arange(0.25, 5.0, 0.02)
for r in radiusRange:
    #storing results to each list
    VResults.append(V(0, r))
    VZeroResults.append(Vzero(r))

#comparison of the results to ensure accuracy and check for fractional error
print('my method', np.sum(VResults))
print('V with z=0', np.sum(VZeroResults))

#plotting the 2 methods to check if it is accurate
fig1, overlay = plt.subplots(figsize=(8, 6))
overlay.plot(radiusRange, VZeroResults, label='eq(4) results')
overlay.plot(radiusRange, VResults, dashes=[6,2], label='eq(3) results')
overlay.set(xlabel='radius (mm)', ylabel='Voltage (V)',
       title='Overlay of V(r, z=0) and V(r,0)')
overlay.grid()
overlay.legend()
fig1.savefig("q4.png", dpi=150)
plt.show()

#part b
#creating a range of z from -5mm to 5mm with 200 entries
zRange = np.arange(-5, 5, 0.05)
#creating a 2D array that is currently filled with 0 with zRange rows and radiusRange columns
VzResults = [[0 for z in range(len(zRange))] for r in range(len(radiusRange))]

for r in range(len(radiusRange)):
    for z in range(len(zRange)):
        #update the 2D array respective values
        VzResults[r][z] = V(zRange[z], radiusRange[r])

#plot contour graph
X, Y = np.meshgrid(zRange, radiusRange)

fig2, contour = plt.subplots(figsize=(10, 8))
CS = contour.contour(X, Y, VzResults, extent=(-5, 5, -0, 5))
contour.set(xlabel='z range (mm)', ylabel='radius range (mm)',
       title='Contour graph of Voltage (V)')
CB = fig2.colorbar(CS, shrink=0.8, extend='both')
contour.clabel(CS, inline=1, fontsize=10)
contour.grid()
fig2.savefig("q4_2.png", dpi=150)
plt.show()
