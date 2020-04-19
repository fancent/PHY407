"""
The purpose of this code is to find out blah blah blah...
"""
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from gaussxw import gaussxwab

def I(x,temp):
    """
    This function takes blah blah blah
    """
    power = (constants.h * constants.c * (1-x)) / (constants.k * temp)
    return (1-x)**3 / (np.exp(power) - 1)
    
def Efficiency(l1, l2, T, N):
    """
    This function takes blah blah blah
    """
    #transform 
    a = 1 - 1/l1
    b = 1 - 1/l2
    #Gaussian Quadrature
    x,w = gaussxwab(N, a, b)
    results = 0
    for i in range(N):
        results += w[i]*I(x[i],T)
    return results

def findEfficiency(l1, l2, T, N):
    """
    This function takes blah blah blah
    """
    return Efficiency(l1,l2,T,N)/Efficiency(1e-8, 1000000, T, N)

#1 data point for every K in the given range
temperatureRange = np.arange(300, 10000, 1)
#results for light-bulb
results = findEfficiency(380e-9,780e-9,temperatureRange, 400)*100
#find the maximum efficiency to within 1K
maximum = np.where(results == np.max(results))
print("maximum efficiency is:", np.max(results), "%", "with temperature", ' '.join(map(str, temperatureRange[maximum[0]])), "K")

#plot for part a
fig1, graph = plt.subplots()
graph.plot(temperatureRange, results)
graph.set(xlabel="Temperature (K)", ylabel="Efficiency (%)",
       title="Energy Efficiency of incandescent light bulb \nover range of temperature in visible")
graph.grid()
fig1.savefig("q2_A.png")
plt.show()


#part b
resultsInfraRed = findEfficiency(780e-9,2250e-9,temperatureRange, 400)*100
#find the maximum efficiency to within 1K
maximum = np.where(resultsInfraRed == np.max(resultsInfraRed))
print("maximum efficiency is:", np.max(resultsInfraRed), "%", "with temperature", ' '.join(map(str, temperatureRange[maximum[0]])), "K")

#plot for part b
fig2, graph = plt.subplots()
graph.plot(temperatureRange, resultsInfraRed)
graph.set(xlabel="Temperature (K)", ylabel="Efficiency (%)",
       title="Energy Efficiency of incandescent light bulb \nover range of temperature in infra-red")
graph.grid()
fig2.savefig("q2_C.png")
plt.show()