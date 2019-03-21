# ppiclF
A parallel particle-in-cell library in Fortran.

## Capabilities
* Integration for the system:
           
           d y
           ---  =  ydot, 
           d t
           
           $ \sum_{\forall i}{x_i^{2}} $
           
  where y and ydot are vectors that are entirely
  user defined.
       
* Open MPI parallelization allows millions of equations
  to be solved.
       
* Load balances equations based on spatial position of
  particles.
       
* Allows simple user input of external overlapping mesh
  for interactions between particles and external mesh,
  including interpolation and projection.
       
* Includes optional fast binned parallel nearest neighbor
  search between particles within a user defined distance.
