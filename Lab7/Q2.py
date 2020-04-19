"""
The purpose of this program is to describe the energy levels
for the first 3 states of a quantum harmonic oscillator, as well
as an anharmonic oscillator given particular potential functions
as well as initial conditions regarding the particle
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# Constants
e = constants.e
hbar = constants.hbar
m = constants.m_e
V0 = 50 * e  # initial energy (J)
initial = -10 ** -10 # initial position (m)
last = 10 ** -10 # last position (m)
psi_initial = 0.0 #initial wavefunction
a = 10 ** -11  # Angstrom
N = 1000  # total steps
h = (last - initial) / N

def V(x):
    """
    This function calculates our potential functions
    one given for the harmonic oscillator and one for the anharmonic
    """
    return V0 * x ** 2 / a ** 2
    #return V0 * x ** 4 / a ** 4

def f(r, x, E):
    """
    This next function calculates the first order differentials
    we set our derivative of psi to equal phi, and our derivative of 
    phi to follow from the schrodinger equation
    """
    psi = r[0]
    phi = r[1]
    return np.array([phi, (2 * m / hbar ** 2) * (V(x) - E) * psi], float)


def rk4(E):
    """
    This function uses rK4 method of solving our first order differentials
    """
    r = np.array([psi_initial, 1.0] ,float)
    wave = []
    for x in np.arange(initial, last, h):
        wave.append(r[0])
        k1 = h * f(r, x, E)
        k2 = h * f(r + 0.5 * k1, x + 0.5 * h, E)
        k3 = h * f(r + 0.5 * k2, x + 0.5 * h, E)
        k4 = h * f(r + k3, x + h, E)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return np.array(wave, float)

def q2(E1, E2):
    """
    This function normalizes our wavefunction
    and as such include a way of solving the integral
    """
    target_accuracy = e / 1000
    wave = rk4(E1)
    psi2 = wave[N - 1]
    while abs(E1 - E2) > target_accuracy:
        wave = rk4(E2)
        psi1, psi2 = psi2, wave[N - 1]
        E1, E2 = E2, E2 - psi2 * (E2 - E1) / (psi2 - psi1)
    mod_squared = wave * wave
    integral = h / 3 *(mod_squared[0] + mod_squared[N//2 - 1] + \
            4 * np.sum(mod_squared[1 : N//2 : 2]) + 2 * np.sum(mod_squared[0 : N//2 + 1 : 2]) )

    return E2 / e, wave / np.sqrt(2*integral)


# PART A
# harmonic oscillator
# note: need to uncomment the other return line in V(x)
# =============================================================================
# E0, psi0 = q2(0, 0.5*e)
# E1, psi1 = q2(200*e, 400*e)
# E2, psi2 = q2(500*e, 700*e)
# print("E_0 = {} eV".format(E0))
# print("E_1 = {} eV".format(E1))
# print("E_2 = {} eV".format(E2))
# =============================================================================

# PART B
# anharmonic oscillator
E0, psi0 = q2(0, 0.5*e)
E1, psi1 = q2(400*e, 600*e)
E2, psi2 = q2(900*e, 1100*e)
print("E_0 = {} eV".format(E0))
print("E_1 = {} eV".format(E1))
print("E_2 = {} eV".format(E2))

# PART C
xpoints = np.arange(initial, last, h)
x_range = slice(N // 4 ,  3 * N // 4 , 1)
fig, graph = plt.subplots(figsize=(8,4))
graph.plot(xpoints[x_range], psi0[x_range], 'k', label="first")
graph.plot(xpoints[x_range], psi1[x_range], 'b', label="second")
graph.plot(xpoints[x_range], psi2[x_range], 'g', label="third")
#graph.set(xlabel='x (m)', ylabel='psi (1/m)', title='Anharmonic oscilator')
graph.set(xlabel='x (m)', ylabel='psi (1/m)', title='Harmonic oscilator')
graph.grid()
graph.legend()
fig.savefig("q2.png")
plt.show()