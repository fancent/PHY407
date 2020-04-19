"""
The purpose of this code is to find out the general form of ellipse equation
by performing matrix operations. Then proceed to find perihelion, aphelion, 
period, eccentricity, semi major and minor axes.
"""

import numpy as np
import matplotlib.pyplot as plt

def semiaxis(a,b,c,d,f,g):
    """
    This function takes in the coefficients of the expanded ellispe function,
    and calculate both the semi-major and semi-minor axes in AU.
    These equations are taken from: http://mathworld.wolfram.com/Ellipse.html
    """
    numerator = 2*((a*(f**2))+(c*(d**2))+(g*(b**2))-(2*b*d*f)-(a*c*g))
    aDenominator = ((b**2)-a*c)*(np.sqrt((a-c)**2 + 4*(b**2)) - (a+c))
    bDenominator = ((b**2)-a*c)*(-np.sqrt((a-c)**2 + 4*(b**2)) - (a+c))
    major = np.sqrt(numerator/aDenominator)
    minor = np.sqrt(numerator/bDenominator)
    return major,minor

def period(a,b,l1,v1):
    """
    This function takes semi-major, semi-minor axes, perihelion, and aphelion
    distance to calculate the period of the orbit in years.
    This equation is taken from Lab4 handout.
    """
    return (2*np.pi*a*b)/(l1*v1)

def eccentricity(a, b):
    """
    This function takes semi-major and semi-minor axes to calculate 
    the eccentricity of the orbit.
    These equations are taken from: http://mathworld.wolfram.com/Ellipse.html
    """
    return np.sqrt(1-((b**2)/(a**2)))

def distances(a, b):
    """
    This function takes semi-major and semi-minor axes to calculate the 
    perihelion and aphelion distances in AU.
    This equation is originally taken from Lab4 handout, since we know a and b,
    we can solve the system of equations to find l1 and l2.
    """
    l1 = a - np.sqrt(a**2 - b**2)
    l2 = a + np.sqrt(a**2 - b**2)
    return l1,l2

def maxVelocity(l1, l2, v2):
    """
    This function takes perihelion, aphelion distances and minimum velocity
    to calculate the maximum velocity at the aphelion in AU/year.
    This equation is taken from Lab4 handout.
    """
    return (l2*v2)/l1


# X and Y given in table 1 in AU
x = [-38.04, -35.28, -25.58, -28.80, -30.06]
y = [27.71, 17.27, 30.68, 31.50, 10.53]

#Calculating parameters for our ellipse equation
#Matrix for X parameters
XMatrix = np.asarray([[x[i]**2, x[i]*y[i], y[i]**2, x[i], y[i], 1] for i in range(len(x))])
#Matrix for Y
YMatrix = np.zeros((6,6))
YMatrix[2][0] = 2
YMatrix[0][2] = 2
YMatrix[1][1] = -1/4
#find the coefficients by calculating the eigenvalues of the following
S_inv = np.linalg.inv(np.dot(XMatrix.T, XMatrix))
E, V = np.linalg.eig(np.dot(S_inv,YMatrix))
maxarg = np.argmax(np.abs(E))
#our ellipse coefficients are inside this list. 
#Note that B, D, F are all multiplied by 2
A,B,C,D,F,G = V[:, maxarg]

#calculating answers to the question
minSpeed = 1350.* (31536000/149598073000) #m/s -> AU/year
semiMajorAxis, semiMinorAxis = semiaxis(A,B/2,C,D/2,F/2,G) #AU
L1, L2 = distances(semiMajorAxis, semiMinorAxis) #AU
maxSpeed = maxVelocity(L1, L2, minSpeed) #AU/year
e = eccentricity(semiMajorAxis, semiMinorAxis)
T = period(semiMajorAxis, semiMinorAxis, L1, maxSpeed) #year
print("Perihelion distance=", L1, "AU")
print("Aphelion distance=", L2, "AU")
print("Period =", T, "years")
print("Eccentricity =", e)
print("Semi-major axis =", semiMajorAxis, "AU")
print("Semi-minor axis =", semiMinorAxis, "AU")

#setting up for the contour graph
#setting the range for our x and y coordinates
xlinspace = np.linspace(-40,5, 400)
ylinspace = np.linspace(-5,35, 400)
xmesh, ymesh = np.meshgrid(xlinspace, ylinspace)
#using our coefficients from above calculations
#we can use the below equations to fit our data points into the orbit
#the contour function will create a range of z values, however we only want
#z = 0. Thus we will display only the first level which is when z=0
z = A*(xmesh**2) + B*xmesh*ymesh +C*(ymesh**2) + D*xmesh + F*ymesh + G

#creating plots
fig1, graph = plt.subplots()
graph.plot(x, y, "x", label="given data points")
graph.plot(0, 0, "o", label="Sun")
graph.contour(xmesh, ymesh, (z), [0])
graph.set(xlabel="x coordinate (AU)", ylabel="y coordinate (AU)",
       title="Orbit of our comet")
graph.grid()
graph.legend()
fig1.savefig("q1.png")
plt.show()