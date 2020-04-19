"""
This file aims to explore on how to use rk4 algorithm to evaluate a system of
ODEs. We will take the SEIR model as our example which simulates the spread of 
infectious disease through a population. This model divides the population into
the categories Susceptible(S), Exposed(E), Infectious(I), and Recovered (R).
We will later introduce the Vaccinated population(V)
"""

import numpy as np
import matplotlib.pyplot as plt

#setting up the constants provided from the lab manuel
B = 0.2
r = 0.1
b = 10**-4
d = 10**-4
a = 0.1
v = 10**-3
N = 1
populations = np.array([1 - 10**-6, 0., 10**-6, 0.],float)

def SEIRmodel(population, t, infectionDeathRate):
    """
    This function represents the system of ODEs for the SEIR model
    """
    S, E, I, R = population
    N = (S+E+I+R)
    lamb = B*I/N
    dSdt = b*N - lamb*S + v*R - d*S
    dEdt = lamb*S - a*E - d*E
    dIdt = a*E - r*I - infectionDeathRate*I
    dRdt = r*I - v*R - d*R
    return np.array([dSdt, dEdt, dIdt, dRdt])
    
def rk4(f, initial, a, b, N, infectionDeathRate):  
    """
    This function performs rk4 to solve for the SEIR's system of ODEs.
    It takes the initial population, start date, end date, number of slices,
    and the mortality rate for the infected population
    """
    h = (b-a)/N
    time = np.arange(a,b,h)
    S_results = []
    E_results = []
    I_results = []
    R_results = []
    for t in time:
        S_results.append(initial[0])
        E_results.append(initial[1])
        I_results.append(initial[2])
        R_results.append(initial[3])
        k1 = h*f(initial, t, infectionDeathRate)
        k2 = h*f(initial+(0.5*k1), t + 0.5*h, infectionDeathRate)
        k3 = h*f(initial+(0.5*k2), t + 0.5*h, infectionDeathRate)
        k4 = h*f(initial + k3, t + h, infectionDeathRate)
        initial += (k1 + (2*k2) + (2*k3) + k4)/6
    return np.asarray(S_results), np.asarray(E_results), np.asarray(I_results), np.asarray(R_results), time

# PART C
# Birth rate = death rate
S, E, I, R, T = rk4(SEIRmodel, populations, 0, 3650, 1000, d)
total = np.asarray(S) + np.asarray(E) + np.asarray(I) + np.asarray(R)
#plotting
fig, graph = plt.subplots()
graph.plot(T, S, label="S")
graph.plot(T, E, label="E")
graph.plot(T, I, label="I")
graph.plot(T, R, label="R")
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title='Population groups over time')
graph.grid()
graph.legend()
fig.savefig("q2c_groups.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(T, total)
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title='Total Population over time')
graph.grid()
fig.savefig("q2c_total.png")
plt.show()


# PART D
# death rate for infected is 100 times more
S, E, I, R, T = rk4(SEIRmodel, populations, 0, 36500, 5000, 10**-2)
total = np.asarray(S) + np.asarray(E) + np.asarray(I) + np.asarray(R)
#plotting
fig, graph = plt.subplots()
graph.plot(T, S, label="S")
graph.plot(T, E, label="E")
graph.plot(T, I, label="I")
graph.plot(T, R, label="R")
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title="Population groups over time with\nlarger infection mortality rate")
graph.grid()
graph.legend()
fig.savefig("q2d_groups.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(T, total, label="Total population")
graph.hlines(0.5, 0, 36500, colors='red', label="50%")
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title="Total Population over time with\nlarger infection mortality rate")
graph.grid()
graph.legend()
fig.savefig("q2d_total.png")
plt.show()

# PART E
# adding a vaccinated population into consideration
vacinnatedPopulations = np.array([1 - 10**-6, 0., 10**-6, 0., 0.],float)
def P(t):
    """birth rate of vaccinated babies"""
    return 1 - 2**(-(t-5475)/365)

def SEIRmodelNoVaccine(population,t, infectionDeathRate):
    """
    same SEIR model as above but added the Vaccinated population(V)
    which we kept its change as constant since the vaccine has not be invented
    """
    S, E, I, R, V = population
    N = (S+E+I+R)
    lamb = B*I/N
    dSdt = b*N - lamb*S + v*R - d*S
    dEdt = lamb*S - a*E - d*E
    dIdt = a*E - r*I - infectionDeathRate*I
    dRdt = r*I - v*R - d*R
    dVdt = 0
    return np.array([dSdt, dEdt, dIdt, dRdt, dVdt])

def SEIRmodelVaccine(population,t, infectionDeathRate):
    """
    same SEIR model as above but added the Vaccinated population(dVdt)
    into considerations, which the new born are vaccinated.
    """
    S, E, I, R, V = population
    N = (S+E+I+R+V)
    p = P(t)
    lamb = B*I/N
    dVdt = p*b*N - d*V
    dSdt = (1-p)*b*N - lamb*S + v*R - d*S
    dEdt = lamb*S - a*E - d*E
    dIdt = a*E - r*I - infectionDeathRate*I
    dRdt = r*I - v*R - d*R
    return np.array([dSdt, dEdt, dIdt, dRdt, dVdt])
    
def rk4_Vaccine(f1, f2, initial, a1, b1, N1, b2, N2, infectionDeathRate):  
    """
    performs rk4 in 2 time interval which the first does not have vaccine
    and the second one has a vaccinated newborn population
    """
    h1 = (b1-a1)/N1
    time15years = np.arange(a1,b1,h1)
    h2 = (b2-b1)/N2
    timeRest = np.arange(b1,b2,h2)
    totalTime = np.concatenate((time15years, timeRest))
    S_results = []
    E_results = []
    I_results = []
    R_results = []
    V_results = []
    
    for t in time15years:
        S_results.append(initial[0])
        E_results.append(initial[1])
        I_results.append(initial[2])
        R_results.append(initial[3])
        V_results.append(initial[4])
        k1 = h1*f1(initial, t, infectionDeathRate)
        k2 = h1*f1(initial+(0.5*k1), t + 0.5*h1, infectionDeathRate)
        k3 = h1*f1(initial+(0.5*k2), t + 0.5*h1, infectionDeathRate)
        k4 = h1*f1(initial + k3, t + h1, infectionDeathRate)
        initial += (k1 + (2*k2) + (2*k3) + k4)/6
        
    for t in timeRest:
        S_results.append(initial[0])
        E_results.append(initial[1])
        I_results.append(initial[2])
        R_results.append(initial[3])
        V_results.append(initial[4])
        k1 = h1*f2(initial, t, infectionDeathRate)
        k2 = h1*f2(initial+(0.5*k1), t + 0.5*h1, infectionDeathRate)
        k3 = h1*f2(initial+(0.5*k2), t + 0.5*h1, infectionDeathRate)
        k4 = h1*f2(initial + k3, t + h1, infectionDeathRate)
        initial += (k1 + (2*k2) + (2*k3) + k4)/6
    return S_results, E_results, I_results, R_results, V_results, totalTime

S, E, I, R, V, T = rk4_Vaccine(SEIRmodelNoVaccine, SEIRmodelVaccine, vacinnatedPopulations, 0, 5475, 5000, 36500, 10000, 10**-2)
total = np.asarray(S) + np.asarray(E) + np.asarray(I) + np.asarray(R) + np.asarray(V)

# plotting
fig, graph = plt.subplots()
graph.plot(T, S, label="S")
graph.plot(T, E, label="E")
graph.plot(T, I, label="I")
graph.plot(T, R, label="R")
graph.plot(T, V, label="V")
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title="Population groups over time with larger infection\nmortality rate and vaccinated population")
graph.grid()
graph.legend()
fig.savefig("q2e_groups.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(T, total, label="Total population")
graph.set(xlabel='time (day)', ylabel='Population (people)',
       title="Total Population over time with larger infection\nmortality rate and vaccinated population")
graph.grid()

fig.savefig("q2e_total.png")
plt.show()

