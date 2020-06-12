# Monte Carlo Solutions to PDE
This project uses the idea of random processes to estimate results of multidimensional partial differential equations, in particular, this project focuses on 2D Laplace's Equation.

I picked a **Monte Carlo approach** to break down the necessary steps of estimating the results asynchronisly.

The details of the implementation is actually very simple as documented in the python notebook file.

The results were an astonishing **99% decrease in time** with an **average error of 3%**.

The results comparison can be seen below (Left: normal computation, Right: Monte Carlo method):

<img src="https://github.com/fancent/PHY407/blob/master/Project/Lab8Method3D.png" width="50%"><img src="https://github.com/fancent/PHY407/blob/master/Project/MonteMethod3D.png" width="50%">

At the end of my project, I briefly explained the approach to tackle PDE in higher complexity and higher dimensions. For full details, please refer to my pdf report.
