import numpy as np
import random
import matplotlib.pyplot as plt

def totalEnergy(arr,J):
    rows = arr[:-1,:]*arr[1:,:]
    cols = arr[:,:-1]*arr[:,1:]
    energy = -J*(np.sum(rows)+np.sum(cols))
    return energy

def sim(N, arr, J, T):
    M = [np.sum(arr)]
    size = arr.shape[0]
    for _ in range(N):
        old_arr = arr.copy()
        old_energy = totalEnergy(old_arr, J)
        rand_i, rand_j = random.randrange(size),random.randrange(size)
        arr[rand_i][rand_j] = - arr[rand_i][rand_j]
        new_energy = totalEnergy(arr, J)
        if new_energy > old_energy:
            if random.random() >= np.exp((old_energy-new_energy)/T):
                arr = old_arr
        M.append(np.sum(arr))
    return M



values = [-1, 1]
size = 20
N = 1000000
temp = 3
a = random.choices(values, k=size**2)
a = np.array(a).reshape(size,size)
results = sim(N, a, 1, temp)

fig, graph = plt.subplots()
graph.plot(results, 'x')
graph.set(xlabel='time step', ylabel='total magnetization',
       title="Total magnetization over {} steps at Temperature={}K".format(N, temp))
graph.grid()
fig.savefig("q2c_{}.png".format(temp))
plt.show()