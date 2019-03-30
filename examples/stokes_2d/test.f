!-----------------------------------------------------------------------
c main code below
!-----------------------------------------------------------------------
      program test
#include "PPICLF"
      include 'mpif.h' 

      call MPI_INIT(ierr) 
      icomm = MPI_COMM_WORLD
      call MPI_COMM_RANK(ppiclf_comm, nid, ierr) 
      call MPI_COMM_SIZE(ppiclf_comm, np , ierr)
 
      call ppiclf_comm_InitMPI(icomm,nid,np)
         call PlaceParticle(npart,ppiclf_y)
      call ppiclf_solve_InitParticle(1,2,0,npart,ppiclf_y) 

      ! time loop
      iostep = 1E2
      nstep  = 1E4
      dt     = 1E-4
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

      integer   npart
      real      y(*)
      real      ran2
      external  ran2

      npart   = 250    ! particles/rank to distribute
      rdum    = ran2(-1-ppiclf_nid) ! initialize random number generator

      do i=1,npart
         ! set initial conditions for solution
         j = PPICLF_LRS*(i-1)
         y(PPICLF_JX +j) = ran2(2)
         y(PPICLF_JY +j) = ran2(2)
         y(PPICLF_JVX+j) = 0.0
         y(PPICLF_JVY+j) = 0.0
      
         ! set some initial particle properties
         ppiclf_rprop(PPICLF_R_JTAUP,i) = 1./9.8
      enddo

      return
      end
!-----------------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
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
