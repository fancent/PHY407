import numpy as np
import matplotlib.pyplot as plt

data = np.unpackbits(np.load("Earth.npy")).reshape(2160, 4320)

#constants
minTheta = 0
maxTheta = 2*np.pi
minPhi = 0
maxPhi = 2*np.pi

#random theta and phi generator
def gen(num):
    phi_result = []
    theta_result = []
    for _ in range(num):
        phi = np.random.uniform(0,2*np.pi)
        phi_result.append(phi)
        theta = np.arcsin(2*np.random.random()-1) + np.pi/2
        theta_result.append(theta)
    return np.array(theta_result), np.array(phi_result)

# thetas, phis = gen(500)
# plt.plot(np.degrees(phis),np.degrees(thetas),'x')
# plt.show()

def land_frac():
    land_result = 0
    water_result = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] == 1:
                land_result += data[i][j] * (np.sin(i*(np.pi/2160))*(np.pi/(2160*4320)))
            else:
                water_result +=  (np.sin(i*(np.pi/2160))*(np.pi/(2160*4320)))
    return land_result/(land_result+water_result)

def land_frac_from_indices(theta_indices, phi_indices):
    land_result = 0
    water_result = 0
    for i in range(len(theta_indices)):
        if data[theta_indices[i]][phi_indices[i]] == 1:
            land_result += (np.sin(i*(np.pi/2160))*(np.pi/(2160*4320)))
        else:
            water_result += (np.sin(i*(np.pi/2160))*(np.pi/(2160*4320)))
    return land_result/(land_result+water_result)

def randomlyCalculateLandFraction(N = 50000):
    # To save computation time I directly set the value from land_frac
    # answer = land_frac()
    answer = 0.2672613705423998
    accuracy = 0.5
    size = 2
    thetas, phis = gen(N)
    istenfound = True
    ispointonefound = True
    iszeropointonefound = True
    while accuracy >= 0.00001 and size <= N:
        size += 10
        theta_index = ((np.pi - thetas) * (2160/np.pi)).astype(int)
        phi_index = (phis * (4320/(2*np.pi))).astype(int)
        test = land_frac_from_indices(theta_index[:size],phi_index[:size])
        accuracy = np.abs(test-answer)/answer
        if istenfound and (accuracy-0.1 <= 0.00001):
            print("land fraction with N={}".format(size),test)
            print("accuracy with N={} at 10%".format(size),accuracy)
            istenfound = False
        if ispointonefound and (accuracy-0.001 <= 0.00001):
            print("land fraction with N={}".format(size),test)
            print("accuracy with N={} at 0.1%".format(size),accuracy)
            ispointonefound = False
        if iszeropointonefound and (accuracy-0.00001 <= 0.00001):
            print("land fraction with N={}".format(size),test)
            print("accuracy with N={} at 0.001%".format(size),accuracy)
            iszeropointonefound = False
    if iszeropointonefound:
        print("accuracy to 0.001% not found")

randomlyCalculateLandFraction()