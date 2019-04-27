[![Build Status](https://travis-ci.org/dpzwick/ppiclF.svg?branch=master)](https://travis-ci.org/dpzwick/ppiclF)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/dpzwick/ppiclF/blob/master/LICENSE)

# ppiclF
ppiclF is a parallel particle-in-cell library written in Fortran. Its main purpose is to provide a unified and scalable interface for a user to solve the following system of differential equations

           
           d y
           ---  =  ydot.
           d t


* See [documentation website](https://dpzwick.github.io/ppiclF-doc) for more details, theory, examples, questions, etc.

* On-the-fly load-balancing of the system of equations across MPI processing ranks based on the coordinates associated with each particle. 

* Simple user input of an external overlapping mesh for interactions between particles and their nearby cells.

* Optional fast binned parallel nearest neighbor search between particles within a user specified distance so that more sophisticated user-implemented right-hand-side forcing models can easily be evaluated. 

* Algorithms have demonstrated scalability to 100,000 processors, allowing billions of equations to be solved simultaneously. 

* Links to both Fortran and C++ external code as a library.
