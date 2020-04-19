import numpy as np
from gaussxw import gaussxw
from scipy import constants
import matplotlib.pyplot as plt

#setting up constants
m1 = 1 #mass (kg)
x0_1 = 0.01 #initial displacement (m)
N0 = 8 #number of data points
N1 = 16 #number of data points
k_s = 12 #spring constant (N/m)
a1 = 0.0 #lower bound
b1 = x0_1 #upper bound

def g(x, x0, k, m):
    halfKX = 0.5 * k * (x0 ** 2 - x ** 2)
    mc2 = m * constants.c ** 2
    numerator = halfKX * (2 * mc2 + halfKX)
    denominator = (mc2 + halfKX)
    return denominator / (constants.c * np.sqrt(numerator))

def T(a,b,N,f,k,x0,m):
    #function that can take different a,b,N to calculate for T
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    results = []
    for i in range(N):
        results.append(4*wp[i]*f(xp[i], x0, k, m))
    return results, f(xp, x0, k, m), xp

integrand8, gResults8, xValues8= T(a1,b1,N0,g,k_s,x0_1,m1)
integrand16, gResults16, xValues16= T(a1,b1,N1,g,k_s,x0_1,m1)
print('sum of integrand for N=8:', np.sum(integrand8))
print('sum of integrand for N=16:', np.sum(integrand16))
print('ideal answer:', 2*np.pi*np.sqrt(m1/k_s))
print("fractional error for N=8", np.abs(np.sum(integrand8) - 2*np.pi*np.sqrt(m1/k_s)))
print("fractional error for N=16", np.abs(np.sum(integrand16) - 2*np.pi*np.sqrt(m1/k_s)))

#Plot for part a N=8
fig1, graph = plt.subplots()
graph.plot(xValues8, gResults8,marker='o',ls='')
graph.set(xlabel="Sample Points from Gaussian quadrature x_k (m)", ylabel="Integrand g_k (m/s)",
       title="Sample Points x_k vs Integrand g_k for N={}".format(N0))
graph.grid()
fig1.savefig("q2_A_1_N={}.png".format(N0))

fig2, graph = plt.subplots()
graph.plot(xValues8, integrand8)
graph.set(xlabel="Sample Points from Gaussian quadrature x_k (m)", ylabel="Weighted Values w_k g_k (s)",
       title="Sample Points x_k vs Weighted Values w_k g_k for N={}".format(N0))
graph.grid()
fig2.savefig("q2_A_2_N={}.png".format(N0))

#Plot for part a N=16
fig3, graph = plt.subplots()
graph.plot(xValues16, gResults16,marker='o',ls='')
graph.set(xlabel="Sample Points from Gaussian quadrature x_k (m)", ylabel="Integrand g_k (m/s)",
       title="Sample Points x_k vs Integrand g_k for N={}".format(N1))
graph.grid()
fig3.savefig("q2_A_1_N={}.png".format(N1))

fig4, graph = plt.subplots()
graph.plot(xValues16, integrand16)
graph.set(xlabel="Sample Points from Gaussian quadrature x_k (m)", ylabel="Weighted Values w_k g_k (s)",
       title="Sample Points x_k vs Weighted Values w_k g_k for N={}".format(N1))
graph.grid()
fig4.savefig("q2_A_2_N={}.png".format(N1))


#part b answer
x0_2 = constants.c/np.sqrt(12) #initial displacement (m)
N2 = 200 #number of data points
a2 = 0.0 #lower bound
b2 = 10*x0_2 #upper bound
#part c
integrand200, gResults200, xValues200 = T(a2,b2,N2,g,k_s,x0_2,m1)
#Plot for part c
fig5, graph = plt.subplots()
graph.plot(xValues200, integrand200)
graph.set(xlabel="Sample Points from Gaussian quadrature x_k (m)", ylabel="Weighted Values w_k g_k (s)",
       title="Sample Points x_k vs Weighted Values w_k g_k for N={}".format(N2))
graph.grid()
fig5.savefig("q2_C_N={}.png".format(N2))
plt.show()