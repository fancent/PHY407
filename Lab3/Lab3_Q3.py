import numpy as np
import matplotlib.pyplot as plt

def myf(x):
    return np.exp(2*x)

def myf3derivative(x):
    return 6*np.exp(2*x)

def delta(f,x,m,h):
    if m> 1:
        return (delta(f, x+h/2, m-1, h) - delta(f, x-h/2, m-1, h))/(h)
    else:
        return (f(x+(h/2))-f(x-(h/2)))/h
    
def centralDiffErr(f,x,C,h):
    return ((2*C*np.abs(f(x)))/h) + (1/24)*(h**2)*(np.abs(delta(myf,x,3,h)))

def optimalError(f,f3,x,C):
    return ((9/8)*(C**2)*(f(x)**2)*f3(x))**(1/3)

#creating a list of number that increase by a factor of 10
hRange = [10**i for i in range(-16,1)]

#using delta to find the first derivative approximation at x = 0
#using central difference method given on the lab 3 document.
#Then find the absolute difference from known solution, 2.
diffFrom2 = [np.abs(delta(myf, 0, 1, i) - 2) for i in hRange]

#graphing part a
fig1, graph = plt.subplots(figsize=(6, 4), dpi=80)
graph.plot(np.arange(-16,1), diffFrom2)
graph.set(xlabel="Power of h (step size)", ylabel="Absolute difference from 2",
       title="Absolute error of each \nderivative over step size")
graph.grid()
fig1.savefig("q3_A.png")

#constant taken from lab manuel
C1 = 10**-16
#calculations for part (a)
print("centtal difference error calculation from textbook",[centralDiffErr(myf, 0, C1, i) for i in hRange])
print("value of h that produced the least difference from 2:", hRange[diffFrom2.index(np.min(diffFrom2))])
print("Min difference from our calculation:", np.min(diffFrom2))
print("Optimal error from equation 5.101:", optimalError(myf, myf3derivative, 0, C1))


#Part B
def cauchy(f,z,m,N):
    coefficients = (float(np.math.factorial(m)))/N
    summation = 0
    for k in range(0, N):
        zk = np.exp(1j*2*np.pi*k/N)
        summation += f(zk)*np.exp(-1j*2*np.pi*k*m/N)
    return (coefficients*summation).real

#Setting up constants and ranges
N = 10000 #number of summation for Cauchy formula
correctDerivative = np.asarray([2**m for m in range(1,11)]) #list of correct result of m-th derivatives
centralDiffResults = np.abs(np.asarray([delta(myf,0,m,0.1) for m in range(1,11)]) - correctDerivative)
cauchyResults = np.abs(np.asarray([cauchy(myf, 0, m, N) for m in range(1,11)]) - correctDerivative)

#graphing part b
fig2, graph = plt.subplots(figsize=(6, 4), dpi=80)
graph.plot(np.arange(1,11), centralDiffResults)
graph.set(xlabel="mth derivative (step size)", ylabel="Absolute difference",
       title="Absolute error of each derivative using\n Central Difference over step size")
graph.grid()
fig2.savefig("q3_B_1.png")

fig3, graph = plt.subplots(figsize=(6, 4), dpi=80)
graph.plot(np.arange(1,11), cauchyResults)
graph.set(xlabel="mth derivative (step size)", ylabel="absolute difference",
       title="Absolute error of each derivative using\n Cauchy Formula over step size")
graph.grid()
fig3.savefig("q3_B_2.png")
