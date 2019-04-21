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

# Statement of Need

Particle-in-cell simulations are ubiquitous in science and engineering. Relevant applications range from explosively driven compressible multiphase gas-solid flow to atomistic modeling. At the core of these diverse applications, a system of differential equations is advanced in time. As a result, the system given by

$$\dfrac{d \mathbf{Y}}{dt} = \dot{\mathbf{Y}}$$

must be solved subject to appropriate initial conditions $\mathbf{Y} (t = 0)$ and right-hand-side forcing $\dot{\mathbf{Y}}$. However, computational constraints make this a daunting task when there are billions of equations, whose forcing may even depend on the states of other equations. Adding to this difficulty is the fact that the particles are often added as an afterthought when building particle-in-cell software. As a result, typical particle-in-cell implementations may initially perform below expectations.

# Summary

ppiclF [@ppiclF] is a parallel particle-in-cell library written in Fortran. Its main purpose is to provide a unified and scalable interface for a user to solve the above general system of differential equations. The library performs on-the-fly load-balancing of the given system of equations across MPI processing ranks based on the coordinates associated with each particle. The library allows simple user input of an external overlapping mesh for interactions between particles and their nearby cells. Additionally, ppiclF includes an optional fast binned parallel nearest neighbor search between particles within a user specified distance so that more sophisticated user-implemented right-hand-side forcing models can easily be evaluated. The algorithms have demonstrated scalability to 100,000 processors, allowing billions of equations to be solved simultaneously. A more complete description of ppiclF and its use is found on the documentation site [@ppiclf-doc].

ppiclF was originally designed to be used by researchers in the field of multiphase flow. In fact, it has been used as a discrete element method solver coupled to both the incompressible Navier-Stokes solver Nek5000 [@nek5000] and the compressible Euler equation solver CMT-nek [@cmt-nek] to solve large-scale fluid-particle systems. Examples of these applications include fluidized beds and multiphase shock-tubes, which are detailed in the dissertation of [@Zwick]. However, since the above system of differential equations are also found in many other disciplines, ppiclF's applications are not limited to multiphase flow. As a result, this library is organized such that a user is only required to set the initial conditions and forcing for the system; other details, such as the parallelization strategy, are handled internally. Thus, ppiclF's framework allows a user to easily translate their mathematical system of equations into executable code.

# Acknowledgements

This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant No. DGE-1315138. This work was also supported in part by the US Department of Energy, National Nuclear Security Administration, Advanced Simulation and Computing Program, as a Cooperative Agreement under the Predictive Science Academic Alliance Program, under Contract No. DE-NA0002378. 

There have been numerous suggestions and reworking of ppiclF over the past few years, including valuable input from the Nek5000 team at Argonne National Laboratory and the Center for Compressible Multiphase Turbulence team at the University of Florida.

# References

