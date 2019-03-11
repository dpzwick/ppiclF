!-----------------------------------------------------------------------
c main code below
!-----------------------------------------------------------------------
      program test
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h' 

      real rparam(lpm_nparam) 
 
      call MPI_INIT(ierr) 
      lpm_comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(lpm_comm, lpm_nid, ierr) 
      call MPI_COMM_SIZE(lpm_comm, lpm_np , ierr)

      rparam(1)  = 1           ! use custom values
      rparam(2)  = 1           ! time integration method
      rparam(3)  = 2           ! polynomial order of mesh
      rparam(4)  = 0           ! use 1 for tracers only
      rparam(5)  = 0.1        ! filter width in real units
      rparam(6)  = 4           ! how many grid points to resolve filter over
      rparam(7)  = 1E-2        ! percent decay of Gaussian filter
      rparam(8)  = 1           ! periodic in x (== 0) ! dont do periodic without bounds!!!
      rparam(9)  = 1           ! periodic in y (== 0)
      rparam(10) = 1           ! periodic in z (== 0)
      rparam(11) = 8E-4        ! time step
      rparam(12) = 3           ! problem dimensions

      call init_particles(lpm_y,npart)
c     call lpm_io_vtu_read('new99999.vtu',npart)
      call lpm_init      (rparam,lpm_y,npart,0.0) 

      ! time loop
      iostep = 100
      nstep  = 1000
      do lpm_cycle=1,nstep
         lpm_time = (lpm_cycle-1)*lpm_dt
         call lpm_solve(lpm_time,lpm_y,lpm_ydot)


         if(mod(lpm_cycle,iostep) .eq. 0) then
             call lpm_io_vtu_write('',0)
c            call lpm_io_vtu_write_bins('',0)
             call lpm_io_vtu_write_grd('',0)
             nptmax = iglmax(lpm_npart,1)
             nptmin = iglmin(lpm_npart,1)
             nptsum = iglsum(lpm_npart,1)
             if (lpm_nid .eq. 0) then
                write(6,'(A,I6,A,E16.10)')  'STEP: ',lpm_cycle,
     >                                    ', TIME: ',lpm_time
                write(6,'(A,I6,A,I6,A,I6)') 'NMAX: ',nptmax,
     >                                    ', NMIN: ',nptmin,
     >                                    ', NAVG: ',nptsum/lpm_np
              endif
         endif
      enddo

      call MPI_FINALIZE(ierr) 

      end
!-----------------------------------------------------------------------
      subroutine init_particles(y,npart)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real      y(*)
      real      pi
      real      ran2
      external  ran2

      npart   = 50       ! particles/rank to distribute
      dp      = 0.01   ! particle diameter
      rhop    = 3307.327 ! particle density
      rdum    = ran2(-1-lpm_nid) ! initialize random number generator
      PI      = 4.D0*DATAN(1.D0)

      do i=1,npart
         ! set initial conditions for solution
         j = LPM_LRS*(i-1)
         y(LPM_JX +j) = 0.1 + 0.8*ran2(2)
         y(LPM_JY +j) = 0.7 + 0.2*ran2(2)
         y(LPM_JZ +j) = 0.1 + 0.8*ran2(2)
         y(LPM_JVX+j) = 0.0
         y(LPM_JVY+j) = 0.0
         y(LPM_JVZ+j) = 0.0
      
         ! set some initial particle properties
         lpm_rprop(LPM_R_JRHOP,i) = rhop
         lpm_rprop(LPM_R_JDP  ,i) = dp
         lpm_rprop(LPM_R_JVOLP,i) = pi/6.0*lpm_rprop(LPM_R_JDP,i)**3
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

