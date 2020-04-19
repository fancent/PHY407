import numpy as np
import matplotlib.pyplot as plt

#Part a
#Text Exercise 6.10 a
#we can modify the number of loops to acheive different levels of accuracy
#in our case we achieved 10^-6 by setting the number of loops to 15
xa=1.0
numOfLoops = 15
for k in range(numOfLoops):
    xa = 1 - np.exp(-2*xa)
    
print("value of x with accuracy 10**-6:", xa)
print("value of 1 - exp(-2*x) with accuracy 10**-6:", 1 - np.exp(-2*xa))
    
#Text Exercise 6.10 b
#similar as above, we first reset the variables and instead of a constant c
#we generated a range of values from 0 to 300 with a stepsize of 0.01
x=1.0
numOfLoops = 15
c = np.arange(0, 3, 0.01)
for k in range(numOfLoops):
    x = 1 - np.exp(-c*x)

fig, graph = plt.subplots()
graph.plot(c, x)
graph.set(xlabel="c", ylabel="x",
       title="Solutions of x for a range of c=(0, 3, 0.01) (Q3a)")
graph.grid()
fig.savefig("q3_A.png")
plt.show()


#Part b
#Text Exercise 6.11 b
def F(c):
    """
    This function uses the similar method as above and find the iterations 
    when it approaches the convergence.
    """
    sigdig = 10 ** -6
    iterations = 1
    def f(x):
        return 1 - np.exp(-c*x)

    def error(x1, x2):
        return (x1 - x2) / (1 - 1 / (c * np.exp(-c * x1)))

    x1 = 1.0 # starting value
    x2 = f(x1)
    while(abs(error(x1, x2)) > sigdig):
        x1, x2 = x2, f(x2)
        iterations += 1
    print('The minimum number of iterations for an accuracy of 10**-6 = ', iterations)
    print("value of x:", x2)
    print("value of 1 - exp(-2*x):", 1 - np.exp(-2*x2))

#Text Exercise 6.11 c
def F_over(c, w):
    """
    This function uses the over-relaxation method to increases the rate of
    convergence by artifically increasing the gradient of f(x)
    """
    sigdig = 10 ** -6
    iterations = 1
    def f(x):
        return 1 - np.exp(-c*x)

    def derivf(x):
        return c * np.exp(-c*x)

    def error(x1, x2):
        return (x1 - x2)/(1 - 1/((1 + w)*derivf(x1) - w))

    x1 = 1.0  # starting value
    x2 = (1 + w) * f(x1) - w * x1
    while abs(error(x1, x2)) > sigdig:
        x1, x2 = x2, (1 + w) * f(x2) - w * x2
        iterations += 1
    print('When we set omega to be ', w,'the minimum number of iterations for an accuracy of 10**-6 is ', iterations)
    print("value of x:", x2)
    print("value of 1 - exp(-2*x):", 1 - np.exp(-2*x2))

print()
print("part b answer")
F(2)
print()
F_over(2, 0.5)


#Part c
sigfig = 10 ** -6
a = 1
b = 2
def f(x, y):
    return y*(a+x**2)

def g(x, y):
    return b/(a+x**2)

def F(x,y):
    return np.sqrt(b/y-a)

def G(x,y):
    return x/(a+x**2)

def solution(f, g):
    iterations = 1
    def delta(x1, x2):
        return (x1 - x2) / x1
    x1 = 0.5
    y1 = 0.5
    x2 = f(x1, y1)
    y2 = g(x1, y1)

    while abs(delta(x1, x2)) > sigfig and abs(delta(y1, y2)) > sigfig:
        if iterations > 1000:
            return 'error'
        
        x1, x2, y1, y2 = x2, f(x2, y2), y2, g(x2, y2)
        iterations += 1
    print('The number of iterations = ', iterations, 'the solution for x ~', x2, 'the solution for y ~', y2)
    return [ x2, y2 ]

print()
print("part c answer")
print(solution(F,G))
