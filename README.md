[![Build Status](https://travis-ci.org/dpzwick/ppiclF.svg?branch=master)](https://travis-ci.org/dpzwick/ppiclF)
[![License: MIT](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/dpzwick/ppiclF/blob/master/LICENSE)
[![status](http://joss.theoj.org/papers/af8a524242ba0c072ef4c42f816b76df/status.svg)](http://joss.theoj.org/papers/af8a524242ba0c072ef4c42f816b76df)

# ppiclF
ppiclF is a parallel particle-in-cell library written in Fortran. 

Applications of ppilcF include element-based particle-in-cell simulations, such as Euler-Lagrange mutliphase flow simulation, immersed boundary methods, and even atomistic-scale modeling. At its essence, ppiclF's main purpose is to provide a unified and scalable interface for a user to solve the following system of differential equations

           
           d y
           ---  =  ydot,
           d t

which are found in all of the previously given particle-in-cell applications. See [documentation website](https://dpzwick.github.io/ppiclF-doc) for more details, theory, examples, questions, etc.

# Capabilities

* On-the-fly load-balancing of the system of equations across MPI processing ranks based on the coordinates associated with each particle. 

* Simple user input of an external overlapping mesh for interactions between particles and their nearby cells.

* Optional fast binned parallel nearest neighbor search between particles within a user specified distance so that more sophisticated user-implemented right-hand-side forcing models can easily be evaluated. 

* Algorithms have demonstrated scalability to 100,000 processors, allowing billions of equations to be solved simultaneously. 

* Links to both Fortran and C++ external code as a library.

# Recommended Citation

| Zwick, D. (2019). ppiclF: A Parallel Particle-In-Cell Library in Fortran. *Journal of Open Source Software.* 4(37), 1400. https://doi.org/10.21105/joss.01400 |
| --- |
