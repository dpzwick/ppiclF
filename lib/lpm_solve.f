!-----------------------------------------------------------------------
      subroutine lpm_init(rparam,y,npart,time_)
#include "LPM"

      real     y(*)
      real     rparam(*)

      lpm_d2chk(2) = 0.0
      lpm_npart = npart

      call lpm_comm_setup

      call lpm_rparam_set(rparam)
      call lpm_tag_init
      call lpm_tag_set
      call lpm_init_filter
      call lpm_domain_size(y,npart)

c     ! get domain bounds
c     call domain_size( lpm_xdrange(1,1),lpm_xdrange(2,1)
c    >                 ,lpm_xdrange(1,2),lpm_xdrange(2,2)
c    >                 ,lpm_xdrange(1,3),lpm_xdrange(2,3))

c     ! send particles to correct rank
      call lpm_interpolate_setup

c     ! two-way coupling init
c     if (int(lpm_rparam(4)) .ne. 1) then
c        call lpm_project
c     endif

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_domain_size(y,npart)
#include "LPM"
      include 'mpif.h'

      real y(*), work(1)

      ! assumed that at least x,y,z in 1,2,3
      lpm_jx = 1
      lpm_jy = 2
      lpm_jz = 3

      xmin =  1E8
      xmax = -1E8
      ymin =  1E8
      ymax = -1E8
      zmin =  1E8
      zmax = -1E8

      do j=1,npart
         i = LPM_LRS*(j-1)
         if (y(lpm_jx+i) .gt. xmax) xmax = y(lpm_jx+i)
         if (y(lpm_jx+i) .lt. xmin) xmin = y(lpm_jx+i)
         if (y(lpm_jy+i) .gt. ymax) ymax = y(lpm_jy+i)
         if (y(lpm_jy+i) .lt. ymin) ymin = y(lpm_jy+i)
         if (y(lpm_jz+i) .gt. zmax) zmax = y(lpm_jz+i)
         if (y(lpm_jz+i) .lt. zmin) zmin = y(lpm_jz+i)
      enddo

      n = 1
      call mpi_allreduce (xmin,work,n,mpi_double_precision,
     >                    mpi_min,lpm_comm,ierr)
      call mpi_allreduce (ymin,work,n,mpi_double_precision,
     >                    mpi_min,lpm_comm,ierr)
      call mpi_allreduce (zmin,work,n,mpi_double_precision,
     >                    mpi_min,lpm_comm,ierr)
      call mpi_allreduce (xmax,work,n,mpi_double_precision,
     >                    mpi_max,lpm_comm,ierr)
      call mpi_allreduce (ymax,work,n,mpi_double_precision,
     >                    mpi_max,lpm_comm,ierr)
      call mpi_allreduce (zmax,work,n,mpi_double_precision,
     >                    mpi_max,lpm_comm,ierr)

      lpm_xdrange(1,1) = xmin
      lpm_xdrange(2,1) = xmin
      lpm_xdrange(1,2) = ymin
      lpm_xdrange(2,2) = ymin
      lpm_xdrange(1,3) = zmin
      lpm_xdrange(2,3) = zmin

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_rparam_set(rparam)
#include "LPM"

      real rparam(*)

c     call rzero(lpm_rparam, lpm_nparam)
      do i=1,lpm_nparam
         lpm_rparam(i) = 0.0
      enddo

      ! set defaults
      if (rparam(1) .eq. 0) then
         lpm_rparam(1)  = 0        ! use custom values
         lpm_rparam(2)  = 1        ! time integration method
         lpm_rparam(3)  = LPM_LX1  ! polynomial order of mesh
         lpm_rparam(4)  = 1        ! use 1 for tracers only
         lpm_rparam(5)  = 0        ! index of filter non-dimensionalization in rprop
         lpm_rparam(6)  = 0        ! non-dimensional Gaussian filter width
         lpm_rparam(7)  = 0        ! percent decay of Gaussian filter
         lpm_rparam(8)  = 1        ! periodic in x
         lpm_rparam(9)  = 1        ! periodic in y
         lpm_rparam(10) = 1        ! periodic in z
         lpm_rparam(11) = -1       ! time step
         lpm_rparam(12) = 3        ! problem dimension

      ! custom values
      else
         do i=1,lpm_nparam
            lpm_rparam(i) = rparam(i)
         enddo

         tol = 1.e-12
         if (wdsize.eq.4) tol = 1.e-6

         if (abs(lpm_rparam(2)) .lt. tol  ) lpm_rparam(2) = 1
         if (abs(lpm_rparam(3)) .lt. tol  ) lpm_rparam(3) = LPM_LX1

      endif

      lpm_dt = lpm_rparam(11)

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_init_filter
#include "LPM"

      rsig  = 0.0
      jdp   = int(lpm_rparam(5))
      filt  = lpm_rparam(6) ! assume already set
      alph  = lpm_rparam(7)

      rsig  = filt*lpm_rprop(jdp,i)/(2.*sqrt(2.*log(2.)))
      lpm_d2chk(2) = rsig*sqrt(-2*log(alph))

c     rdx_max = 0.0

c     lx1m1 = max(1,LPM_LX1-1)
c     ly1m1 = max(1,LPM_LY1-1)
c     lz1m1 = max(1,LPM_LZ1-1)

c     ! first, compute what is the smallest filter size for the grid
c     do ie=1,nelt
c     do k=1,lz1m1
c     do j=1,ly1m1
c     do i=1,lx1m1
c        ip1 = i+1
c        jp1 = j+1
c        kp1 = k+1

c        rdx_test = xm1(ip1,j,k,ie) - xm1(i,j,k,ie)
c        rdy_test = ym1(i,jp1,k,ie) - ym1(i,j,k,ie)
c        if (if3d) 
c    >   rdz_test = zm1(i,j,kp1,ie) - zm1(i,j,k,ie)

c        rdist = max(rdx_test,rdy_test)
c        if (if3d)
c    >   rdist = max(rdz_test,rdist)

c        rdx_dy_test = sqrt(rdx_test**2 + rdy_test**2)
c        rdist = max(rdx_dy_test,rdist)
c        if (if3d) then
c           rdy_dz_test = sqrt(rdy_test**2 + rdz_test**2)
c           rdist = max(rdy_dz_test,rdist)
c           rdx_dz_test = sqrt(rdx_test**2 + rdz_test**2)
c           rdist = max(rdx_dz_test,rdist)
c        endif

c        if (rdist .gt. rdx_max) rdx_max = rdist
c     enddo
c     enddo
c     enddo
c     enddo
c     rdx_max = glmax(rdx_max,1)

c     lpm_rdx_max = rdx_max

c     rsig_dx_min_set  = 0.25
c     rfilt_dp_min_set = 1.0
c     rsig_dx_min      = 1E12

c     tol = 1.e-12
c     if (wdsize.eq.4) tol = 1.e-6

c     rfilt = lpm_rparam(6)
c     lpm_d2chk(2) = -1
c     if (abs(alph) .lt. tol) then
c        alph = 1E-3
c        lpm_rparam(7) = alph
c     endif

c     do i=1,lpm_npart

c        if (filt .lt. tol) then
c           rsig_test_grid = rsig_dx_min_set*rdx_max + tol
c           rsig_test_diam = rfilt_dp_min_set*lpm_rprop(jdp,i)
c    >                       /(2.*sqrt(2.*log(2.)))
c           rsig_test = max(rsig_test_grid,rsig_test_diam)
c           filt      = rsig_test/lpm_rprop(jdp,i)*2.*sqrt(2.*log(2.))
c           rfilt = filt
c        else
c           rsig_test = filt*lpm_rprop(jdp,i)/(2.*sqrt(2.*log(2.)))
c        endif

c        rsig_dx = rsig_test/rdx_max
c        if (rsig_dx .lt. rsig_dx_min) rsig_dx_min = rsig_dx

c        if (rsig_test .gt. rsig) rsig = rsig_test
c     enddo
c     rsig = glmax(rsig,1)
c     rsig_dx_min = glmin(rsig_dx_min,1)
c        
c     lpm_d2chk(2) = rsig*sqrt(-2*log(alph))

c     if (rsig_dx_min .lt. rsig_dx_min_set) then
c        if (nid .eq. 0) then
c           filt_new = filt/(rsig_dx_min/rsig_dx_min_set)
c           write(6,100) filt, filt_new
c100        format('Reset Gaussian filter width from', E14.7
c    >             ' to', E14.7, ' or larger')
c           call exitt
c        endif
c     endif

c     lpm_rparam(6) = glmax(rfilt,1)
c     lpm_d2chk(2)  = glmax(lpm_d2chk(2),1)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_tag_set
#include "LPM"

      do i=1,lpm_npart
         if (lpm_iprop(5,i) .eq. -1) lpm_iprop(5,i) = 0 ! nid
         if (lpm_iprop(6,i) .eq. -1) lpm_iprop(6,i) = 0 ! istep
         if (lpm_iprop(7,i) .eq. -1) lpm_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_tag_init
#include "LPM"

      if (.not. LPM_RESTART) then
         do i=1,LPM_LPART
            lpm_iprop(5,i) = -1
            lpm_iprop(6,i) = -1
            lpm_iprop(7,i) = -1
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_solve(time_,y,ydot)
#include "LPM"

      real time_
      real y(*)
      real ydot(*)

      if (int(lpm_rparam(2)) .eq. 1) call lpm_rk3_driver(time_,y,ydot)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_rk3_driver(time_,y,ydot)
#include "LPM"

      real time_
      real y(*)
      real ydot(*)

      ndum = LPM_NPART*LPM_LRS

      ! save stage 1 solution
      do i=1,ndum
         lpm_y1(i) = y(i)
      enddo

      ! get rk3 coeffs
      call lpm_rk3_coeff

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call lpm_fun(time_,y,ydot)

         ndum = LPM_NPART*LPM_LRS
         ! rk3 integrate
         do i=1,ndum
            y(i) =  rk3coef(1,istage)*lpm_y1 (i)
     >            + rk3coef(2,istage)*y      (i)
     >            + rk3coef(3,istage)*ydot   (i)
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_interpolate_setup
#include "LPM"

c     call lpm_move_outlier
      call lpm_comm_bin_setup
      call lpm_comm_findpts
      call lpm_comm_crystal

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_interpolate_fld(jp,fld)
#include "LPM"

      common /intp_h/ ih_intp(2,1)

      real fld(*)

      ih_intp1 = ih_intp(1,i_fp_hndl)

c     call fgslib_findpts_eval_local( ih_intp1
c    >                               ,lpm_rprop (jp,1)
c    >                               ,LPM_LRP
c    >                               ,lpm_iprop (2,1)
c    >                               ,LPM_LIP
c    >                               ,lpm_rprop2(1,1)
c    >                               ,LPM_LRP2
c    >                               ,LPM_NPART
c    >                               ,fld)

      ! this one
c        call fgslib_findpts_eval( ih_intp1
c    >                           ,lpm_rprop (jp,1)
c    >                           ,LPM_LRP
c    >                           ,lpm_iprop (1,1)
c    >                           ,LPM_LIP
c    >                           ,lpm_iprop (3,1)
c    >                           ,LPM_LIP
c    >                           ,lpm_iprop (2,1)
c    >                           ,LPM_LIP
c    >                           ,lpm_rprop2(1,1)
c    >                           ,LPM_LRP2
c    >                           ,LPM_NPART
c    >                           ,fld)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_project
#include "LPM"

      if (int(lpm_rparam(4)) .ne. 1) then
c        call lpm_comm_ghost_create
c        call lpm_comm_ghost_send
c        call lpm_solve_project
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_rk3_coeff
#include "LPM"

      rk3coef(1,1) = 0.0
      rk3coef(2,1) = 1.0 
      rk3coef(3,1) = lpm_dt
      rk3coef(1,2) = 3.0/4.0
      rk3coef(2,2) = 1.0/4.0 
      rk3coef(3,2) = lpm_dt/4.0 
      rk3coef(1,3) = 1.0/3.0
      rk3coef(2,3) = 2.0/3.0 
      rk3coef(3,3) = lpm_dt*2.0/3.0 

      return
      end
!-----------------------------------------------------------------------
