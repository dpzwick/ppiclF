#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!----------------------------------------------------------------------
! Note that include 'PPICLF' is not allowed except with multiple defs
! linker option in make file, which is not supported on newer mac 
! linkers (-m option, or -z muldefs with gcc)
!----------------------------------------------------------------------
      program main
      include 'mpif.h' 
!
      integer*4 np, nid, icomm, ierr
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

      ! Test A
      call testA
      call test_diagnostic('A')
      call test_diagnostic_timeint('A')

      ! Test B
      call testB
      call test_diagnostic('B')

      ! Test C
      call testC
      call test_diagnostic('C')

      ! Test D
      call testD
      call test_diagnostic('D')

      ! Test E
      call testE
      call test_diagnostic('E')
      call test_diagnostic_projection('E')

      ! Finalize MPI
      call MPI_FINALIZE(ierr) 

      end program
!----------------------------------------------------------------------
      subroutine testA
      include "PPICLF"
!
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) 
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) 
      real*8 time, dt
!
      imethod = 1
      ndim    = 2
      iendian = 0
      npart   = 0
      if (ppiclf_nid .eq. 0) then
         npart = 1

         y(PPICLF_JX,1)  = 0.0d0
         y(PPICLF_JY,1)  = 0.0d0
         y(PPICLF_JVY,1) = 0.0d0
         rprop(PPICLF_R_JTAUP,1) = 1.0d0/9.8d0
      endif
      call ppiclf_solve_InitParam(imethod,ndim,iendian)
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))
!     call ppiclf_solve_InitNeighborBin(1.5d0, 0)

      dt = 1E-4
      do istep=1,1000
         time = (istep-1)*dt
         call ppiclf_solve_IntegrateParticle(istep,istep+1,dt,time)
      enddo

      return
      end
!----------------------------------------------------------------------
      subroutine testB
      include "PPICLF"
!
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) 
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) 
!
      imethod = 1
      ndim    = 2
      iendian = 0
      npart   = 0
      if (ppiclf_nid .eq. 0) then
         npart = 2

         y(PPICLF_JX,1)  = 0.0d0
         y(PPICLF_JY,1)  = 0.0d0
         rprop(PPICLF_R_JTAUP,1) = 1.0d0/9.8d0

         y(PPICLF_JX,2)  = 1.0d0
         y(PPICLF_JY,2)  = 0.0d0
         rprop(PPICLF_R_JTAUP,2) = 1.0d0/9.8d0
      endif
      call ppiclf_solve_InitParam(imethod,ndim,iendian)
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))
!     call ppiclf_solve_InitNeighborBin(1.5d0, 0)

      call ppiclf_solve_IntegrateParticle(1,2,0.0d0,0.0d0)

      return
      end
!----------------------------------------------------------------------
      subroutine testC
      include "PPICLF"
!
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) 
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) 
!
      imethod = 1
      ndim    = 2
      iendian = 0
      npart   = 0
      if (ppiclf_nid .eq. 0) then
         npart = 2

         y(PPICLF_JX,1)  = 0.0d0
         y(PPICLF_JY,1)  = 0.0d0
         rprop(PPICLF_R_JTAUP,1) = 1.0d0/9.8d0

         y(PPICLF_JX,2)  = 1.0d0
         y(PPICLF_JY,2)  = 1.0d0
         rprop(PPICLF_R_JTAUP,2) = 1.0d0/9.8d0
      endif
      call ppiclf_solve_InitParam(imethod,ndim,iendian)
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))
      call ppiclf_solve_InitNeighborBin(1.5d0, 0)

      call ppiclf_solve_IntegrateParticle(1,2,0.0d0,0.0d0)

      return
      end
!----------------------------------------------------------------------
      subroutine testD
      include "PPICLF"
!
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) 
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) 
!
      imethod = 1
      ndim    = 3
      iendian = 0
      npart   = 0
      if (ppiclf_nid .eq. 0) then
         npart = 2

         y(PPICLF_JX,1)  = 0.0d0
         y(PPICLF_JY,1)  = 0.0d0
         y(PPICLF_JZ,1)  = 0.0d0
         rprop(PPICLF_R_JTAUP,1) = 1.0d0/9.8d0

         y(PPICLF_JX,2)  = 1.0d0
         y(PPICLF_JY,2)  = 0.0d0
         y(PPICLF_JZ,2)  = 1.0d0
         rprop(PPICLF_R_JTAUP,2) = 1.0d0/9.8d0
      endif
      call ppiclf_solve_InitParam(imethod,ndim,iendian)
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))
      call ppiclf_solve_InitNeighborBin(1.5d0, 0)

      call ppiclf_solve_IntegrateParticle(1,2,0.0d0,0.0d0)

      return
      end
!----------------------------------------------------------------------
      subroutine testE
      include "PPICLF"
!
      integer*4 imethod, ndim, iendian, npart
      real*8 y(PPICLF_LRS    , PPICLF_LPART) 
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) 
!
      imethod = 1
      ndim    = 3
      iendian = 0
      npart   = 0
      if (ppiclf_nid .eq. 0) then
         npart = 3

         y(PPICLF_JX,1)  = 0.0d0
         y(PPICLF_JY,1)  = 0.0d0
         y(PPICLF_JZ,1)  = 0.0d0
         rprop(PPICLF_R_JTAUP,1) = 1.0d0/9.8d0

         y(PPICLF_JX,2)  = 1.0d0
         y(PPICLF_JY,2)  = 0.0d0
         y(PPICLF_JZ,2)  = 1.0d0
         rprop(PPICLF_R_JTAUP,2) = 1.0d0/9.8d0

         y(PPICLF_JX,3)  = 1.0d0
         y(PPICLF_JY,3)  = 1.0d0
         y(PPICLF_JZ,3)  = 1.0d0
         rprop(PPICLF_R_JTAUP,3) = 1.0d0/9.8d0
      endif
      call ppiclf_solve_InitParam(imethod,ndim,iendian)
      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               npart     ,
     >                               y(1,1)    ,
     >                               rprop(1,1))
      call ppiclf_solve_InitNeighborBin(1.5d0, 0)
      call ppiclf_solve_InitGaussianFilter(0.25d0,1E-2, 0)

      call ppiclf_solve_IntegrateParticle(1,2,0.0d0,0.0d0)


      return
      end
!----------------------------------------------------------------------
      subroutine test_diagnostic_timeint(casename)
      include "PPICLF"
!
      character*1  casename
      character*4  caseflag
      real*8 rint, rerr
      real*8 ppiclf_glmin
      external ppiclf_glmin
!
      rint = ppiclf_glmin(ppiclf_y(PPICLF_JVY,1),1)
      rerr = abs(rint - (dexp(-0.1d0*9.8d0) - 1.0d0)) ! 1 particle (t = 0.1, taup=g)

      write(caseflag,'(4A1)') '(', casename, ')', ' '

      if (ppiclf_nid .eq. 0) then
         write(6,'(A4,A14,E14.8)') caseflag, 'TimeI error:  ', rerr
      endif

      return
      end
!----------------------------------------------------------------------
      subroutine test_diagnostic_projection(casename)
      include "PPICLF"
!
      character*1  casename
      character*4  caseflag
      real*8 rint, rcc, rvol, rerr
      real*8 ppiclf_glsum
      external ppiclf_glsum
!
      ppiclf_ngrids = 6
      call ppiclf_comm_CreateSubBin
      call ppiclf_solve_ProjectParticleSubBin
c     call ppiclf_io_WriteSubBinVTU('')

      ! loop over cells
      ifld = PPICLF_P_JTEST
      rint = 0.0d0
      navg = 4
      if (ppiclf_ndim .eq. 3) navg = 8
      nmax = max(1,ppiclf_bz-1)
      do k=1,nmax
      do j=1,ppiclf_by-1
      do i=1,ppiclf_bx-1
         rcc = ppiclf_grid_fld(i  ,j  ,k  ,ifld) +
     >         ppiclf_grid_fld(i+1,j  ,k  ,ifld) +
     >         ppiclf_grid_fld(i  ,j+1,k  ,ifld) +
     >         ppiclf_grid_fld(i+1,j+1,k  ,ifld)
         if (ppiclf_ndim .eq. 3)
     >   rcc = rcc                               +
     >         ppiclf_grid_fld(i  ,j  ,k+1,ifld) +
     >         ppiclf_grid_fld(i+1,j  ,k+1,ifld) +
     >         ppiclf_grid_fld(i  ,j+1,k+1,ifld) +
     >         ppiclf_grid_fld(i+1,j+1,k+1,ifld)
         rcc = rcc / real(navg) ! take cell average 
         rvol = ppiclf_rdx*ppiclf_rdy ! mult by volume cell
         if (ppiclf_ndim .eq. 3)
     >   rvol = rvol*ppiclf_rdz
         rint = rint + rcc*rvol
      enddo
      enddo
      enddo
      rint = ppiclf_glsum(rint,1)
      rerr = abs(rint - 3.0d0/9.8d0) ! 3 particles project 1/taup

      write(caseflag,'(4A1)') '(', casename, ')', ' '

      if (ppiclf_nid .eq. 0) then
         write(6,'(A4,A14,E14.8)') caseflag, 'Prjct error:  ', rerr
      endif

      return
      end
!----------------------------------------------------------------------
      subroutine test_diagnostic(casename)
      include "PPICLF"
!
      character*1  casename
      character*4  caseflag
      integer*4   npart_total , npart_min , npart_max
      integer*4   nghost_total, nghost_min, nghost_max
      integer*4   ppiclf_iglsum, ppiclf_iglmin, ppiclf_iglmax
      external    ppiclf_iglsum, ppiclf_iglmin, ppiclf_iglmax
!
      npart_total  = ppiclf_iglsum(ppiclf_npart   ,1)
      npart_min    = ppiclf_iglmin(ppiclf_npart   ,1)
      npart_max    = ppiclf_iglmax(ppiclf_npart   ,1)
      nghost_total = ppiclf_iglsum(ppiclf_npart_gp,1)
      nghost_min   = ppiclf_iglmin(ppiclf_npart_gp,1)
      nghost_max   = ppiclf_iglmax(ppiclf_npart_gp,1)

      write(caseflag,'(4A1)') '(', casename, ')', ' '

      if (ppiclf_nid .eq. 0) then
         write(6,'(A4,A14,I2)') caseflag, 'npart total:  ', npart_total
         write(6,'(A4,A14,I2)') caseflag, 'npart min:    ', npart_min
         write(6,'(A4,A14,I2)') caseflag, 'npart max:    ', npart_max
         write(6,'(A4,A14,I2)') caseflag, 'nghost total: ', nghost_total
         write(6,'(A4,A14,I2)') caseflag, 'nghost min:   ', nghost_min
         write(6,'(A4,A14,I2)') caseflag, 'nghost max:   ', nghost_max
         write(6,'(A4,A14,I2)') caseflag, 'nbin x:       ',
     >                                                  ppiclf_n_bins(1)
         write(6,'(A4,A14,I2)') caseflag, 'nbin y:       ',
     >                                                  ppiclf_n_bins(2)
         write(6,'(A4,A14,I2)') caseflag, 'nbin z:       ',
     >                                                  ppiclf_n_bins(3)
         write(6,'(A4,A14,I2)') caseflag, 'nbin total:   ', 
     >                                                  ppiclf_n_bins(1)
     >                                                 *ppiclf_n_bins(2)
     >                                                 *ppiclf_n_bins(3)
      endif

      return
      end
!----------------------------------------------------------------------
