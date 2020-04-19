"""
This file aims to solve the Burger's Equation
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
eps = 1
dx = 0.02
dt = 0.005
L_x = 2*np.pi
T_f = 2.
N_x = int(L_x/dx)
N_t = int(T_f/dt)



# Create arrays
u = np.zeros([N_x+1,N_t+1],float)

# Setting initial boundary conditions
u[:,0] = np.array([np.sin(i*dx) for i in range(len(u))])
u[0] = 0.
u[-1] = 0.

# One forward time step
u[:,1] = u[:,0] - np.array([np.sin(i*dx)*np.cos(i*dx) for i in range(len(u))])*dt
tend = 0.5

# Main loop
t = 0.0
beta = eps/2 * (dt/dx)

for j in range(1,N_t):
    for i in range(N_x):
        if (i != 0 and N_x) and j >= 1:
            u[i][j+1] = u[i][j-1] - beta*(u[i+1][j]**2 - u[i-1][j]**2)
        else:
            u[i][j+1] = u[i][j-1]


fig = plt.figure()
img = plt.imshow(u)
cb = plt.colorbar(img,fraction=0.076, pad=0.04)
cb.set_label('Temperature (C)')
plt.title('Temperature distribution in a heat conductor')
plt.xlabel('x coordinates (cm)')
plt.ylabel('y coordinates (cm)')
plt.savefig("q3.jpg")
plt.show()
plt.figure()
plt.plot(u[:,0])
plt.plot(u[:,100])
plt.plot(u[:,200])
plt.plot(u[:,250])
plt.title('Temperature distribution in a heat conductor')
plt.xlabel('x coordinates (cm)')
plt.ylabel('y coordinates (cm)')
plt.savefig("q3_1.jpg")