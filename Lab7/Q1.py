"""
This file aims to simulate the interaction between 2 Argon
atoms affected by the Lennard-Jones potential with the
Verlet Method
"""
import numpy as np
import matplotlib.pyplot as plt

def r(r1,r2):
    """
    This function calculates the distance between r1 and r2.
    """
    return np.linalg.norm(r1-r2)

def f(r,dx,dy):
    """
    This function calculates the Lennard-Jones potential equations
    """
    force = 48*r**-13 - 24*r**-7
    xforce = force*(dx/r)
    yforce = force*(dy/r)
    if xforce >= 0 or yforce >= 0:
        return np.array([[xforce, yforce], [-xforce, -yforce]])
    else:
        return np.array([[-xforce, -yforce], [xforce, yforce]])

def verlet(f, dt, v0, r1, r2, endtime):
    """
    This function uses the Verlet method to keep track of the 
    positions of r1 and r2 under the Lennard-Jones potential.
    """ 
    x1List = []
    y1List = []
    x2List = []
    y2List = []
    dtRange = np.arange(0, endtime, dt)
    position = np.array([r1,r2])
    r0 = r(r1,r2)
    distDiff = position[0]-position[1]
    v_t_halfh = v0 + 0.5*dt*f(r0, distDiff[0], distDiff[1])
    for _ in dtRange:
        position += dt*v_t_halfh
        dist = r(position[0], position[1])
        distDiff = position[0]-position[1]
        k = dt*f(dist, distDiff[0], distDiff[1])
        #v_t_h = v_t_halfh + 0.5*k we dont need this for this question
        v_t_halfh += k
        x1List.append(position[0][0])
        y1List.append(position[0][1])
        x2List.append(position[1][0])
        y2List.append(position[1][1])
    return x1List, y1List, x2List, y2List, dtRange


a = np.array([4, 4])
b = np.array([5.2, 4])
x1, y1, x2, y2, time = verlet(f, 0.01, 0.0, a, b, 1)

#plotting
fig, graph = plt.subplots()
graph.plot(x1, y1, '.')
graph.plot(x2, y2, '.')
graph.set(xlabel='x coordinate', ylabel='y coordinate',
       title='Trajectories of 2 argon particles')
graph.grid()
fig.savefig("q1_1.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(time, x1, '.')
graph.set(xlabel='time', ylabel='x position',
       title='x position over time for r1')
graph.grid()
fig.savefig("q1_1_xaxis.png")
plt.show()

a = np.array([4.9, 4])
b = np.array([5.1, 4])
x1, y1, x2, y2, time = verlet(f, 0.01, 0.0, a, b, 1)

#plotting
fig, graph = plt.subplots()
graph.plot(x1, y1, '.')
graph.plot(x2, y2, '.')
graph.set(xlabel='x coordinate', ylabel='y coordinate',
       title='Trajectories of 2 argon particles')
graph.grid()
fig.savefig("q1_2.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(time, x1, '.')
graph.set(xlabel='time', ylabel='x position',
       title='x position over time for r1')
graph.grid()
fig.savefig("q1_2_xaxis.png")
plt.show()

a = np.array([2., 3.])
b = np.array([7., 6.])
x1, y1, x2, y2, time = verlet(f, 0.01, 0.0, a, b, 100)

#plotting
fig, graph = plt.subplots()
graph.plot(x1, y1, '.')
graph.plot(x2, y2, '.')
graph.set(xlabel='x coordinate', ylabel='y coordinate',
       title='Trajectories of 2 argon particles')
graph.grid()
fig.savefig("q1_3.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(time, x1, '.')
graph.set(xlabel='time', ylabel='x position',
       title='x position over time for r1')
graph.grid()
fig.savefig("q1_3_xaxis.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(time, y1, '.')
graph.set(xlabel='time', ylabel='y position',
       title='y position over time for r1')
graph.grid()
fig.savefig("q1_3_yaxis.png")
plt.show()