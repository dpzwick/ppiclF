#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
      program main
      include 'mpif.h' 
!
      integer*4 np, nid, icomm
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) ! Normal ordering
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) ! Normal ordering

      integer*4 nstep, iostep
      real*8 dt, time
      integer*4 ierr
      real*8 rdum, ran2
      external ran2
!
      ! Init MPI
      call MPI_Init(ierr) 
      icomm = MPI_COMM_WORLD
      call MPI_Comm_rank(icomm, nid, ierr) 
      call MPI_Comm_size(icomm, np , ierr)

      ! Pass to library to Init MPI
      call ppiclf_comm_InitMPI(icomm,
     >                         nid  ,
     >                         np   )

      ! Set initial conditions and parameters for particles
      imethod = 1
      ndim    = 2
      iendian = 0
      npart   = 250
      rdum    = ran2(-1-nid) ! init random numbers
      do i=1,npart

         y(PPICLF_JX,i)  = ran2(2)
         y(PPICLF_JY,i)  = ran2(2)
         y(PPICLF_JVX,i) = 0.0
         y(PPICLF_JVY,i) = 0.0

         rprop(PPICLF_R_JTAUP,i) = 1.0/9.8
      enddo
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))

      ! Integrate particles in time
      nstep  = 1000
      iostep = 100
      dt     = 1E-4
      do istep=1,nstep

         time = (istep-1)*dt
         call ppiclf_solve_IntegrateParticle(istep ,
     >                                       iostep,
     >                                       dt    ,
     >                                       time  )
      enddo

      ! Finalize MPI
      call MPI_FINALIZE(ierr) 

      stop
      end
!-----------------------------------------------------------------------
      real*8 FUNCTION ran2(idum)
      INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
! Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
! Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
! between 0.0 and 1.0 (exclusive of the endpoint values). 
! Call with idum a negative integer to initialize; thereafter, do not alter 
! idum between successive deviates in a sequence. RNMX should approximate the 
! largest floating value that is less than 1.
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
