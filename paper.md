---
title: 'ppiclF: A Parallel Particle-In-Cell Library in Fortran'
tags:
  - particle
  - particle-in-cell
  - ODE integrator
  - Euler-Lagrange
  - discrete element method
authors:
  - name: David Zwick
    orcid: 0000-0001-5329-6341
    affiliation: 1
affiliations:
 - name: Department of Mechanical & Aerospace Engineering, University of Florida
   index: 1
date: 19 April 2019
bibliography: paper.bib
---

# Summary

Particle-in-cell simulations are ubiquitous in science and engineering. Relevant applications range from explosively driven compressible multiphase gas-solid flow to atomistic modeling. At the core of these diverse applications, a system of differential equations is advanced in time. As a result, the system given by

$$\dfrac{d \mathbf{Y}}{dt} = \dot{\mathbf{Y}}$$

must be solved subject to appropriate initial conditions and right-hand-side forcing. However, computational constraints make this a daunting task when there are billions of equations which may even depend on other equations. Adding to this difficulty is the fact that the particles are often added as an afterthought when building particle-in-cell software. As a result, particle-in-cell implementations may initially perform below expectations.

ppiclF [@ppiclF] is a parallel particle-in-cell library written in Fortran. Its main purpose is to provide a unified and scalable interface for a user to solve the above general system of differential equations. The library performs on-the-fly load-balancing of the given system of equations across MPI processing ranks based on the coordinates associated with each particle. The library allows simple user input of an external overlapping mesh for interactions between particles and nearby cells. Additionally, ppiclF includes an optional fast binned parallel nearest neighbor search between particles within a user defined distance so that more sophisticated user-implemented right-hand-side forcing models can easily be evaluated. The algorithms have been shown to be scalable to 100,000 processors while solving billions of equations. A more complete description of ppiclF and its use is found on the documentation site [@ppiclf-doc].

ppiclF was originally designed to be used by researchers in the field of compressible multiphase flow. It has largely come about as a result of the dissertation of [@Zwick].

# Acknowledgements

This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant No. DGE-1315138. This work was also supported in part by the US Department of Energy, National Nuclear Security Administration, Advanced Simulation and Computing Program, as a Cooperative Agreement under the Predictive Science Academic Alliance Program, under Contract No. DE-NA0002378. 

There have been numerous suggestions and reworking of ppiclF over the past few years, including valuable input from the Nek5000 team at Argonne National Laboratory and the Center for Compressible Multiphase Turbulence team at the University of Florida.

# References

