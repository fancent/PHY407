import numpy as np

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

print(solution(F,G))
    
