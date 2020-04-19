"""
This file aims to simulate the Belousov–Zhabotinsky reaction, 
is a chemical mixture which, when heated, undergoes a series
of reactions that cause the chemical concentrations in the 
mixture to oscillate between two extremes (x, y).
"""
import numpy as np
import matplotlib.pyplot as plt

## Constants 
a = 1
b = 3
x0 = 0  # x concentration level (M)
y0 = 0  # y concentration level (M)
targetAcc = 10 ** -10  # target accuracy for BS method
start = 0 # start time (s)
end = 20 # end time (s)

def f(r):
    """
    This function calculates the equations for the BZ reaction
    """
    x = r[0]
    y = r[1]
    dxdt = 1 - ((b + 1) * x) + (a * (x ** 2) * y)
    dydt = (b * x) - (a * (x ** 2) * y)
    return np.array([dxdt, dydt], float)
    
def midpoint(r, n, H):
    """
    This function calculates the modified mid-point method
    given in the textbook
    """
    r2 = np.copy(r)
    h = H / n
    r1 = r + 0.5 * h * f(r)
    r2 += h * f(r1)
    for _ in range(n - 1):
        r1 += h * f(r2)
        r2 += h * f(r1)
    return 0.5 * (r1 + r2 + 0.5 * h * f(r2))


def BZ_reaction():
    """
    This function simulates the entire Belousov–Zhabotinsky reaction
    from start time to end time with the given constants at the 
    beginning of the file using the Bulirsch–Stoer method with 
    recursion instead of a while loop.
    """
    r = np.array([x0, y0], float)
    tpoints = [start]
    xpoints = [r[0]]
    ypoints = [r[1]]

    def BS(r, t, H):
        """
        This function is just a shell for the following recursive
        function if n, the number of recursive calls, exceeds 8. 
        Then we will redo the calculation with a smaller H. 
        """
        def BS_row(R1, n):
            """
            This function calculates the row of extrapolation estimates.
            Then it calculates the error and check if it falls under 
            our desired accuracy. If not, it will recurse on itself
            with a larger n. If yes, then it will update the list of 
            variables.
            """
            if n > 8:
                r1 = BS(r, t, H / 2)
                return BS(r1, t + H / 2, H / 2)
            else:
                R2 = [midpoint(r, n, H)]
                for m in range(1, n):
                    R2.append(R2[m - 1] + (R2[m - 1] - R1[m - 1]) / ((n / (n - 1)) ** (2 * (m)) - 1))

                R2 = np.array(R2, float)
                error_vector = (R2[n - 2] - R1[n - 2]) / ((n / (n - 1)) ** (2 * (n - 1)) - 1)
                error = np.sqrt(error_vector[0] ** 2 + error_vector[1] ** 2)

                target_accuracy = H * targetAcc
                if error < target_accuracy:
                    tpoints.append(t + H)
                    xpoints.append(R2[n - 1][0])
                    ypoints.append(R2[n - 1][1])
                    return R2[n - 1]
                else:
                    return BS_row(R2, n + 1)     
        return BS_row(np.array([midpoint(r, 1, H)], float), 2)
    
    BS(r, start, end - start)
    return tpoints, xpoints, ypoints

#plotting our results
t, x, y = BZ_reaction()
fig, graph = plt.subplots()
graph.plot(t, x, 'r', label="x")
graph.plot(t, y, 'b', label="y")
graph.plot(t, x, 'r.')
graph.plot(t, y, 'b.')
graph.set(xlabel='time (s)', ylabel='concentration level (M)',
       title='Belousov–Zhabotinsky concentration level over time')
graph.grid()
graph.legend()
fig.savefig("q3.png")
plt.show()