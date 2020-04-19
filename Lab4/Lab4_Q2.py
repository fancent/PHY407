"""
The purpose of this code is to find the maximum efficiency of light bulb
through out a range of temperatures in Kelvin and specific wavelengths.
"""
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from gaussxw import gaussxwab

def zfunc(x):
    """
    This function maps 0 -> inf into 0 -> 1
    """
    return x/(1-x)
    
def I(x,temp):
    """
    This function takes temperature and wavelength and calculate the intensities
    """
    coefficients = 2*np.pi*constants.h*(constants.c**2)
    power = (constants.h * constants.c) / (x*constants.k * temp)
    return coefficients * (x**-5)/(np.exp(power) - 1)

def Ix(x, temp):
    """
    This function maps the above function from 0 -> inf into 0 -> 1
    """
    return I(zfunc(x), temp) * (1/(1-x)**2)

def Efficiency(l1, l2, T, N, f):
    """
    This function uses Gaussian Quadrature method to calculate the integrals
    needed for the efficiency given the bounds of wavelengths.
    """
    x,w = gaussxwab(N, l1, l2)
    results = 0
    for i in range(N):
        results += w[i]*f(x[i],T)
    return results

def findEfficiency(l1, l2, T, N):
    """
    This function calculates the efficiency of the light bulb with
    """
    return Efficiency(l1,l2,T,N,I)/Efficiency(0, 1, T, N,Ix)


#1 data point for every K in the given range
temperatureRange = np.arange(300, 10000, 1)
#results for light-bulb efficiencies
results = findEfficiency(380e-9,780e-9,temperatureRange, 10000)*100

#plotting for part a
fig1, graph = plt.subplots()
graph.plot(temperatureRange, results)
graph.set(xlabel="Temperature (K)", ylabel="Efficiency (%)",
       title="Energy Efficiency of incandescent light bulb \nover range of temperature in visible")
graph.grid()
fig1.savefig("q2_A.png")
plt.show()

#PART B
#find the maximum efficiency to within 1K
maximum = np.where(results == np.max(results))
print("maximum efficiency is:", np.max(results), "%", "with temperature", ' '.join(map(str, temperatureRange[maximum[0]])), "K")

#PART C
#calculates the efficiency of the light bulb in infra-red spectrum
resultsInfraRed = findEfficiency(780e-9,2250e-9,temperatureRange, 10000)*100
#find the maximum efficiency to within 1K
maximum = np.where(resultsInfraRed == np.max(resultsInfraRed))
print("maximum efficiency is:", np.max(resultsInfraRed), "%", "with temperature", ' '.join(map(str, temperatureRange[maximum[0]])), "K")
fig2, graph = plt.subplots()
graph.plot(temperatureRange, resultsInfraRed)
graph.set(xlabel="Temperature (K)", ylabel="Efficiency (%)",
       title="Energy Efficiency of incandescent light bulb \nover range of temperature in infra-red")
graph.grid()
fig2.savefig("q2_C.png")
plt.show()