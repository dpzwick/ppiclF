!-----------------------------------------------------------------------
c main code below
c#include "PPICLF"
!-----------------------------------------------------------------------
      program test
      include 'mpif.h' 

      call MPI_INIT(ierr) 
      ppiclf_comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(ppiclf_comm, ppiclf_nid, ierr) 
      call MPI_COMM_SIZE(ppiclf_comm, ppiclf_np , ierr)

      call ppiclf_comm_InitMPI(ppiclf_comm,ppiclf_nid,ppiclf_np)
c        call PlaceParticle(npart,ppiclf_y)
c     call ppiclf_solve_InitParticle(1,3,0,npart,ppiclf_y) 
c     call ppiclf_solve_InitNeighborBin(0.07)
c     call ppiclf_solve_InitWall(0.0,1.0,0.0,0.0,0.0,0.0)

      ! time loop
      iostep = 100
      nstep  = 1000
      dt     = 1E-2
      do istep=1,nstep
         time = (istep-1)*dt
c        call ppiclf_solve_IntegrateParticle(istep,iostep,dt,time
c    >                                      ,ppiclf_y,ppiclf_ydot)
      enddo

      call MPI_FINALIZE(ierr) 

      end
!-----------------------------------------------------------------------
