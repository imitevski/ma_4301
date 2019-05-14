# Numerical Methods for Monge-Ampere Equation
This work investigates multigrid methods for solving Monge-Ampere equation. We exploit the monotonicity of the equation to write it in an alternative way and then solve it numerically with the Full Approximation Scheme. 

## Abstract: 
The Monge-Ampere (MA) equation is a fully nonlinear degenerate elliptic partial differential equation that arises in optimal mass transportation, beam shaping, image registration, seismology, etc. In the classical form this equation is given by $\det(D^2\phi(x)) = f(x)$ where $\phi$ is constrained to be convex. Previous work has produced solvers that are fast but can fail on realistic (non-smooth) data or robust but relatively slow. The purpose of this work is to implement a more robust and time-efficient scheme for solving the MA equation and do convergence studies for different discretization and full multigrid schemes. We express the MA operator as the product of the eigenvalues of the Hessian matrix. This allows for a globally elliptic discretization that is provably convergent. The method combines a nonlinear Gauss-Seidel iterative method with different discretizations which is stable because the underlying scheme preserves monotonicity. In order to solve these systems efficiently, the V-cycle full approximation scheme multigrid method is exploited with error correction within the recursive algorithm; this scheme is used to leverage the low cost of computation on the coarse grids to build up the finer grids. This work shows computational results that demonstrate the speed and robustness of the algorithm.

## Table of Contents: 
* `code_matlab`: Folder contains the MATLAB version of the code used. Code is run from `main.m`.
* `figures`: Figures used in the report
* `ma_solver.ipynb`: Jupyter Notebook containing the Python version of the code which at this moment has bugs and it will be completed in future
* `report.pdf`: Report submitted for the class for which this project was mentioned 

## Acknowledgements:
The MATLAB code was developed by: [Matthew Illingworth](https://github.com/octoract), Ivan Mitevski, and David Yousuf during EXTREEMS-QED project in 2017 at New Jersey Institute of Technology. Abstract submitted to JMM 2018, San Diego, CA can be found [here](http://jointmathematicsmeetings.org/amsmtgs/2197_abstracts/1135-65-1870.pdf).
