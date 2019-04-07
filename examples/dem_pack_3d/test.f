!-----------------------------------------------------------------------
c main code below
!-----------------------------------------------------------------------
      program test
#include "PPICLF.h"
#include "PPICLF"
      include 'mpif.h' 
!
! Internal:
!
      integer*4 icomm, nid, np
      integer*4 istep, iostep, npart
      real*8 dt, time
      ! For user implemented collision model
      real*8 ksp,erest
      common /external_user_collsion/ ksp,erest
      ! For user implemented collision model
!
      call MPI_INIT(ierr) 
      icomm = MPI_COMM_WORLD
      call MPI_COMM_RANK(icomm, nid, ierr) 
      call MPI_COMM_SIZE(icomm, np , ierr)

      call ppiclf_comm_InitMPI(icomm,nid,np)
         call PlaceParticle(npart,ppiclf_y)
      call ppiclf_solve_InitParticle(1,3,0,npart,ppiclf_y) 
      call ppiclf_solve_InitNeighborBin(0.07)

      call ppiclf_io_ReadWallVTK("geometry/ppiclf_tank.vtk")

      ! For user implemented collision model
      ksp      = 10000.0
      erest    = 0.1
      rpi      = 4.0*atan(1.0)
      rmij1    = rpi/6.0*(0.07**3)*2500. ! change with diff ic's
      rmij2    = rpi/6.0*(0.05**3)*2500. ! change with diff ic's
      nres     = 20
      rmij     = 1./(1./rmij2 + 1./rmij2)
      dt_c_max = sqrt(rmij/ksp*(log(erest)**2 + rpi**2))/nres
      ! For user implemented collision model

      ! time loop
      iostep = 25
      nstep  = 20000
      dt     = dt_c_max
      do istep=1,nstep
         time = (istep-1)*dt
         call ppiclf_solve_IntegrateParticle(istep,iostep,dt,time
     >                                      ,ppiclf_y,ppiclf_ydot)
      enddo

      call MPI_FINALIZE(ierr) 

      stop
      end
!-----------------------------------------------------------------------
      subroutine PlaceParticle(npart,y)
#include "PPICLF"

      integer*4 npart
      real*8    y(*)
      real*8    pi
      real*8    ran2
      external  ran2

      npart   = 200    ! particles/rank to distribute
      dp_min  = 0.05   ! particle diameter min
      dp_max  = 0.07   ! particle diameter max
      rhop    = 2500. ! particle density
      rdum    = ran2(-1-ppiclf_nid) ! initialize random number generator
      PI      = 4.D0*DATAN(1.D0)

      do i=1,npart
         ! set initial conditions for solution
         j = PPICLF_LRS*(i-1)
         y(PPICLF_JX +j) = -0.2  + 0.4*ran2(2)
         y(PPICLF_JY +j) = -0.2  + 0.6*ran2(2)
         y(PPICLF_JZ +j) = -0.2  + 0.4*ran2(2)
         y(PPICLF_JVX+j) = 0.0
         y(PPICLF_JVY+j) = 0.0
         y(PPICLF_JVZ+j) = 0.0
      
         ! set some initial particle properties
         ppiclf_rprop(PPICLF_R_JRHOP,i) = rhop
         ppiclf_rprop(PPICLF_R_JDP  ,i) = dp_min+(dp_max-dp_min)*ran2(2)
         ppiclf_rprop(PPICLF_R_JVOLP,i) = pi/6.0
     >                                  *ppiclf_rprop(PPICLF_R_JDP,i)**3
      enddo

      return
      end
!-----------------------------------------------------------------------
      real*8 FUNCTION ran2(idum)
      INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER*4 idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then 
         idum1=max(-idum,1) 
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1 
            if (idum1.lt.0) idum1=idum1+IM1 
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1) 
      endif
      k=idum1/IQ1 
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1 
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2 
      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1 
      if(iy.lt.1)iy=iy+IMM1 
      ran2=min(AM*iy,RNMX)
      return
      END
c----------------------------------------------------------------------
