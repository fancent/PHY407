"""
This function aims to simulate Temperature distribution in a heat conductor
with a given initial boundary conditions. We will be using Gauss-Seidel 
relaxation with replacement and overrelaxation to calculate the temperatures.
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
M = 200 # width of boundary (mm)
h = 80 # height of boundary (mm)
target = 1e-6   # Target accuracy

# Create arrays to hold temperature values
phi = np.zeros([h+1,M+1],float)

# Setting initial boundary conditions
phi[40,90:110] = 25 #C to D
phi[0,:] = 20 #G to H
for i in range(90):
    phi[80,110+i] = 5 - 5/90 * i # E to F
    phi[80,90-i] = 5 - 5/90 * i # B to A  
for i in range(80):
    phi[0+i,0] = 20 - 20/80 * i # H to A
    phi[0+i,200] = 20 - 20/80 * i # G to F
for i in range(40):
    phi[40+i,90] = 25 - 20/40 * i # H to A
    phi[40+i,110] = 25 - 20/40 * i # G to F

omega = 0.9 # relaxation coefficient
delta = 2 * target # initial accuracy

# The following loop updates each pixel value via Gauss-Seidel 
# relaxation with replacement and overrelaxation.
# Please comment the while loop and uncomment for loop for part (b)

while delta>target:
#for _ in range(100):
    delta = 0.0
    for i in range(h):
        for j in range(M):
            if (i!=0 and i!=h and j!=0 and j!=M) and (not (40<=i<=80 and 90<=j<=110)):
                oldvalue = phi[i, j]
                phi[i, j] = (1 + omega) * (phi[i + 1, j] + phi[i - 1, j] + phi[i, j + 1] + phi[i, j - 1]) / 4 \
                          - omega * phi[i, j]

                delta = max(delta, np.abs(phi[i, j] - oldvalue))
                
print("temperature at x=2.5cm y=1cm:", phi[55][10],"C")

# Make a plot
fig = plt.figure()
img = plt.imshow(phi, extent=[0,20,0,8])
cb = plt.colorbar(img,fraction=0.076, pad=0.04)
cb.set_label('Temperature (C)')
plt.title('Temperature distribution in a heat conductor')
plt.xlabel('x coordinates (cm)')
plt.ylabel('y coordinates (cm)')
plt.savefig("q1d_w={}.jpg".format(omega))
plt.show()
