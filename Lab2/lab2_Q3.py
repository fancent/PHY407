import numpy as np
from scipy import constants

#defining the function to be integrated with
def f(x):
    return (x**3)/((np.exp(x)) - 1.)

#defining the simpsons rule calculation
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

#constants for simpsons rule where N is the number of slices which must be even
# a is the lower bound which we picked a really small number to approximate 0 to 
# avoid reaching division by zero
# b is the upper bound which we picked a big number to approximate infinity 
# since we are going to the e^700 which is reaching python's limit and we do not
# want to have overflow issues
N = 10000
a = 0.000001
b = 700

#checking integral with wolfram alpha's result
integral = simpsonsRule(f, a, b, N)
print('integral:', integral)

#let temperature to be 100 Kelvin for our calculations
T = 100
#The constant that we got from part a
C = (2 * constants.pi * (constants.k**4) * (T**4))/((constants.h**3)*(constants.c**2))
W = C * integral
#comparing our results and checking the accuracy
print('constant from integration', W/(T**4))
print('scipy constant', constants.sigma)
print('Accuracy', (1 - ((W/(T**4)) - constants.sigma)/constants.sigma) * 100, '%')
