[![Build Status](https://travis-ci.org/dpzwick/ppiclF.svg?branch=master)](https://travis-ci.org/dpzwick/ppiclF)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/dpzwick/ppiclF/blob/master/LICENSE)

# ppiclF
A parallel particle-in-cell library in Fortran.

* Integration for the system:
           
           d y
           ---  =  ydot, 
           d t
           
  where y and ydot are vectors that are entirely
  user defined.
       
* Open MPI parallelization allows billions of equations
  to be solved.
       
* Load balances equations based on spatial position of
  particles.

* Links with both Fortran and C++ external code.
       
* Allows simple user input of external overlapping mesh
  for interactions between particles and external mesh,
  including interpolation and projection.
       
* Includes optional fast binned parallel nearest neighbor
  search between particles within a user defined distance.
