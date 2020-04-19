import numpy as np
import matplotlib.pyplot as plt

#define function of p(u)
def p(u):
    return (1-u)**8

#define function of q(u)
def q(u):
    return (1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 -56*u**5 +28*u**6 -8*u**7 +u**8)

#create data set for (a) of about 500 elements
data = np.arange(0.98, 1.02, 0.00008)

#part a
#initialize following lists to store results from p(u), q(u), 
#and p(u) - q(u) for graph plotting
pResult = []
qResult = []
pqDiffResult = []

#going through each data point and adding the results to respective lists
for i in data:
    pResult.append(p(i))
    qResult.append(q(i))
    pqDiffResult.append(p(i) - q(i))

#plotting the results of p(u) and q(u)
fig1, pGraph = plt.subplots()
pGraph.plot(data, pResult)
pGraph.set(xlabel='range of values around 1', ylabel='result of p function',
       title='p function graph')
pGraph.grid()
fig1.savefig("q2_1.png")
plt.show()

fig2, qGraph = plt.subplots()
qGraph.plot(data, qResult)
qGraph.set(xlabel='range of values around 1', ylabel='result of q function',
       title='q function graph')
qGraph.grid()
fig2.savefig("q2_2.png")
plt.show()


#part b
#calculating terms for eq 7
C = 10**-16
#we treat exponents similar to multiplication, N = number of operations
N = 8 + (8 + 7 + 6 + 5 + 4 + 3 + 2 + 1)
#we treat sum of x similar to sum of coefficients since u to very close to 1
coefficients = [1, 8, 28, 56, 70, 56, 28, 8, 1]
sumOfSquared = 0
for i in coefficients:
    sumOfSquared += i**2
    
#using eq 7 to computer standard deviation
std7 = C*np.sqrt(sumOfSquared)
#reference answer using numpy.std method
diffStd = np.std(pqDiffResult, ddof=1)
#compare their values
print('numpy std method:', diffStd)
print('eq 7 method:',std7)

#plotting histogram
n_bins = 18
histogram, axs = plt.subplots()
axs.hist(pqDiffResult, bins=n_bins)
axs.set(xlabel='p(u) - q(u)', ylabel='Frequency',
       title='Histogram of p(u) - q(u)')
histogram.savefig("q2_histogram.png")
plt.show()


#part c
#create data set of 0.980 to 0.984 with about 200 entries
dataHalf = np.arange(0.980, 0.984, 0.00002)
#initialize list for storing the results
pqDiffHalfRangeResult = []
for i in dataHalf:
    #finding values of u where the fractional error is close to 100%
    if 0.98 <= (np.abs(p(i) - q(i))/np.abs(p(i))) <= 1.02:
        print('factional error is about 100% at', i)
    #storing the results to the list
    pqDiffHalfRangeResult.append(np.abs(p(i) - q(i))/np.abs(p(i)))

#plotting the scattered plot
scat, scatter = plt.subplots()
scatter.plot(dataHalf, pqDiffHalfRangeResult, "o", markersize=5, label='fractional error of each value')
scatter.hlines(1.0, 0.980, 0.984, color='red', label='100% error')
scatter.set(xlabel='range', ylabel='result of abs(p-q)/abs(p)',
       title='abs(p-q)/abs(p) from 0.980 to close to 1.0')
scatter.grid()
scatter.legend()
scat.savefig("q2_scat.png")
plt.show()


#part d
#defining function for part d
def f(u):
    return (u**8/((u**4)*(u**4)))
#initialize list for storing the results
fResult = []
for i in data:
    #storing results to the list
    fResult.append(f(i)-1)
    
#plotting
fig3, fGraph = plt.subplots()
fGraph.plot(data, fResult)
fGraph.set(xlabel='data', ylabel='result of f function',
       title='f function graph')
fGraph.grid()
fig3.savefig("q2_3.png")
plt.show()

