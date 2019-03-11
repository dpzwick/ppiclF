!-----------------------------------------------------------------------
      subroutine lpm_init(rparam,y,npart,time_)
#include "lpm_user.h"
#include "lpm.h"
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

c     ! get domain bounds
      lpm_xdrange(1,1) = -1E8
      lpm_xdrange(2,1) =  1E8
      lpm_xdrange(1,2) = -1E8
      lpm_xdrange(2,2) =  1E8
      lpm_xdrange(1,3) = -1E8
      lpm_xdrange(2,3) =  1E8
c     call domain_size( lpm_xdrange(1,1),lpm_xdrange(2,1)
c    >                 ,lpm_xdrange(1,2),lpm_xdrange(2,2)
c    >                 ,lpm_xdrange(1,3),lpm_xdrange(2,3))

      lpm_nbmax = LPM_BMAX ! set only one bin for now...

c     ! send particles to correct rank
      call lpm_interpolate_setup

c     ! two-way coupling init
      call lpm_project

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_rparam_set(rparam)
#include "lpm_user.h"
#include "lpm.h"
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
         lpm_rparam(5)  = 0        ! filter size in real units

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
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      rsig  = 0.0
      alph  = lpm_rparam(7)

      rsig  = lpm_rparam(5)/(2.*sqrt(2.*log(2.)))
      lpm_d2chk(2) = rsig*sqrt(-2*log(alph))

c     if (int(lpm_rparam(4)) .eq. 1) lpm_d2chk(2) = 0

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
#include "lpm_user.h"
#include "lpm.h"
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
#include "lpm_user.h"
#include "lpm.h"
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
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real time_
      real y(*)
      real ydot(*)

      if (int(lpm_rparam(2)) .eq. 1) call lpm_rk3_driver(time_,y,ydot)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_rk3_driver(time_,y,ydot)
#include "lpm_user.h"
#include "lpm.h"
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
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      call lpm_move_outlier
      call lpm_comm_bin_setup
      call lpm_comm_findpts
      call lpm_comm_crystal

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_interpolate_fld(jp,fld)
#include "lpm_user.h"
#include "lpm.h"
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
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      if (int(lpm_rparam(4)) .ne. 1) then
         call lpm_comm_ghost_create
         call lpm_comm_ghost_send
         call lpm_solve_project
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_rk3_coeff
#include "lpm_user.h"
#include "lpm.h"
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
      subroutine lpm_move_outlier
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      integer in_part(LPM_LPART), jj(3), iperiodicx, iperiodicy,
     >                                   iperiodicz

      iperiodicx = int(lpm_rparam(8))
      iperiodicy = int(lpm_rparam(9))
      iperiodicz = int(lpm_rparam(10))
      ndim       = int(lpm_rparam(12))

      jj(1) = 1
      jj(2) = 2
      jj(3) = 3

      do i=1,lpm_npart
         isl = (i -1) * LPM_LRS + 1
         in_part(i) = 0
         do j=0,ndim-1
            jchk = jj(j+1)
            if (lpm_y(jchk,i).lt.lpm_xdrange(1,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   lpm_y(jchk,i) = lpm_xdrange(2,j+1) - 
     &                         abs(lpm_xdrange(1,j+1) - lpm_y(jchk,i))
                   lpm_y1(isl+j)   = lpm_xdrange(2,j+1) +
     &                         abs(lpm_xdrange(1,j+1) - lpm_y1(isl+j))
                  goto 1512
                endif
            endif
            if (lpm_y(jchk,i).gt.lpm_xdrange(2,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   lpm_y(jchk,i) = lpm_xdrange(1,j+1) +
     &                         abs(lpm_y(jchk,i) - lpm_xdrange(2,j+1))
                   lpm_y1(isl+j)   = lpm_xdrange(1,j+1) +
     &                         abs(lpm_y1(isl+j) - lpm_xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (lpm_iprop(1,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1512 continue
         enddo
      enddo

      ic = 0
      do i=1,lpm_npart
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               isl = (i -1) * LPM_LRS + 1
               isr = (ic-1) * LPM_LRS + 1
               call copy(lpm_y     (1,ic),lpm_y(1,i)     ,LPM_LRS)
               call copy(lpm_y1    (isr) ,lpm_y1(isl)    ,LPM_LRS)
               call copy(lpm_ydot  (1,ic),lpm_ydot(1,i)  ,LPM_LRS)
               call copy(lpm_ydotc (1,ic),lpm_ydotc(1,i) ,LPM_LRS)
               call copy(lpm_rprop (1,ic),lpm_rprop(1,i) ,LPM_LRP)
               call copy(lpm_rprop2(1,ic),lpm_rprop2(1,i),LPM_LRP2)
               call icopy(lpm_iprop(1,ic),lpm_iprop(1,i) ,LPM_LIP)
            endif
         endif
      enddo
      lpm_npart = ic

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_solve_project
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real    multfci
      integer e

      real    rproj(1+LPM_LRP_GP,LPM_LPART+LPM_LPART_GP)
      integer iproj(4,LPM_LPART+LPM_LPART_GP)

      logical partl

      real pi

      PI=4.D0*DATAN(1.D0)

      nxyz = LPM_BX1*LPM_BY1*LPM_BZ1

      nxyzdum = nxyz*LPM_LRP_PRO
      do i=1,nxyzdum
         lpm_grid_fld(i,1,1,1) = 0.0
      enddo

      d2chk2_sq = lpm_d2chk(2)**2

      ! real particles
      lpm_jxgp  = 1
      lpm_jygp  = 2
      lpm_jzgp  = 3

c     lpm_npart_gp = 0

      ! real particles
      do ip=1,lpm_npart
         rsig    = lpm_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (lpm_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip) = rbexpi
         rproj(2 ,ip) = lpm_cp_map(lpm_jxgp,ip)
         rproj(3 ,ip) = lpm_cp_map(lpm_jygp,ip)
         rproj(4 ,ip) = lpm_cp_map(lpm_jzgp,ip)


         do j=5,LPM_LRP_GP+1
            rproj(j,ip) = lpm_cp_map(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - lpm_binx(1,1))/lpm_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - lpm_biny(1,1))/lpm_rdy)
         iproj(3,ip) = 
     >       floor( (rproj(4,ip) - lpm_binz(1,1))/lpm_rdz)
      enddo

      ! ghost particles
      do ip=1,lpm_npart_gp
         rsig    = lpm_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (lpm_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip+lpm_npart) = rbexpi
         rproj(2 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jxgp,ip)
         rproj(3 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jygp,ip)
         rproj(4 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jzgp,ip)

         do j=5,LPM_LRP_GP+1
            rproj(j,ip+lpm_npart) = lpm_rprop_gp(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip+lpm_npart) = 
     >       floor( (rproj(2,ip+lpm_npart) - lpm_binx(1,1))/lpm_rdx)
         iproj(2,ip+lpm_npart) = 
     >       floor( (rproj(3,ip+lpm_npart) - lpm_biny(1,1))/lpm_rdy)
         iproj(3,ip+lpm_npart) = 
     >       floor( (rproj(4,ip+lpm_npart) - lpm_binz(1,1))/lpm_rdz)
      enddo

      ndum = lpm_npart+lpm_npart_gp
c     ndum = lpm_npart_gp
c     ndum = lpm_npart

      idum = floor(lpm_rparam(5)/2.0/lpm_rdx
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1
      jdum = floor(lpm_rparam(5)/2.0/lpm_rdy
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1
      kdum = floor(lpm_rparam(5)/2.0/lpm_rdz
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1

      do ip=1,ndum
      ! project here
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_solve_project_bins
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real    multfci
      integer e

      real    rproj(1+LPM_LRP_GP,LPM_LPART+LPM_LPART_GP)
      integer iproj(4,LPM_LPART+LPM_LPART_GP)

      logical partl

      real pi

      PI=4.D0*DATAN(1.D0)

      nxyz = LPM_BX1*LPM_BY1*LPM_BZ1

      nxyzdum = nxyz*LPM_LRP_PRO
      do i=1,nxyzdum
         lpm_grid_fld(i,1,1,1) = 0.0
      enddo

      d2chk2_sq = lpm_d2chk(2)**2

      ! real particles
      lpm_jxgp  = 1
      lpm_jygp  = 2
      lpm_jzgp  = 3

c     lpm_npart_gp = 0

      ! real particles
      do ip=1,lpm_npart
         rsig    = lpm_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (lpm_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip) = rbexpi
         rproj(2 ,ip) = lpm_cp_map(lpm_jxgp,ip)
         rproj(3 ,ip) = lpm_cp_map(lpm_jygp,ip)
         rproj(4 ,ip) = lpm_cp_map(lpm_jzgp,ip)


         do j=5,LPM_LRP_GP+1
            rproj(j,ip) = lpm_cp_map(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - lpm_binx(1,1))/lpm_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - lpm_biny(1,1))/lpm_rdy)
         iproj(3,ip) = 
     >       floor( (rproj(4,ip) - lpm_binz(1,1))/lpm_rdz)
      enddo

      ! ghost particles
      do ip=1,lpm_npart_gp
         rsig    = lpm_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (lpm_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip+lpm_npart) = rbexpi
         rproj(2 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jxgp,ip)
         rproj(3 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jygp,ip)
         rproj(4 ,ip+lpm_npart) = lpm_rprop_gp(lpm_jzgp,ip)

         do j=5,LPM_LRP_GP+1
            rproj(j,ip+lpm_npart) = lpm_rprop_gp(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip+lpm_npart) = 
     >       floor( (rproj(2,ip+lpm_npart) - lpm_binx(1,1))/lpm_rdx)
         iproj(2,ip+lpm_npart) = 
     >       floor( (rproj(3,ip+lpm_npart) - lpm_biny(1,1))/lpm_rdy)
         iproj(3,ip+lpm_npart) = 
     >       floor( (rproj(4,ip+lpm_npart) - lpm_binz(1,1))/lpm_rdz)
      enddo

      ndum = lpm_npart+lpm_npart_gp
c     ndum = lpm_npart_gp
c     ndum = lpm_npart

      idum = floor(lpm_rparam(5)/2.0/lpm_rdx
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1
      jdum = floor(lpm_rparam(5)/2.0/lpm_rdy
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1
      kdum = floor(lpm_rparam(5)/2.0/lpm_rdz
     >    *sqrt(-log(lpm_rparam(7))/log(2.0)))+1

c     if (int(lpm_rparam(3)) .eq. 1) then
         do ip=1,ndum
            iip = iproj(1,ip)
            jjp = iproj(2,ip)
            kkp = iproj(3,ip)

            il  = max(1     ,iip-idum)
            ir  = min(lpm_bx,iip+idum)
            jl  = max(1     ,jjp-jdum)
            jr  = min(lpm_by,jjp+jdum)
            kl  = max(1     ,kkp-kdum)
            kr  = min(lpm_bz,kkp+kdum)

c           write(6,*) il,ir,jl,jr,kl,kr

c           do k=1,lpm_bz
c           do j=1,lpm_by
c           do i=1,lpm_bx
            do k=kl,kr
            do j=jl,jr
            do i=il,ir
    
c              if (lpm_grid_ii(i,j,k) .gt. ir) cycle
c              if (lpm_grid_ii(i,j,k) .lt. il) cycle
c              if (lpm_grid_jj(i,j,k) .gt. jr) cycle
c              if (lpm_grid_jj(i,j,k) .lt. jl) cycle
c              if (lpm_grid_kk(i,j,k) .gt. kr) cycle
c              if (lpm_grid_kk(i,j,k) .lt. kl) cycle


               rdist2  = (lpm_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                   (lpm_grid_y(i,j,k) - rproj(3,ip))**2
               if(lpm_rparam(12) .gt. 2) rdist2 = rdist2 +
     >                   (lpm_grid_z(i,j,k) - rproj(4,ip))**2

               if (rdist2 .gt. d2chk2_sq) cycle

c              write(6,*) 'Made itt', lpm_grid_x(i,j,k)
c    >                              , lpm_grid_y(i,j,k)
c    >                              , lpm_grid_z(i,j,k)

            
               rexp = exp(rdist2*rproj(1,ip))
               
               do jj=1,LPM_LRP_PRO
                  j1 = jj+4
                  lpm_grid_fld(i,j,k,jj) = 
     >                            lpm_grid_fld(i,j,k,jj) 
     >                          + rproj(j1,ip)*rexp
               enddo
            enddo
            enddo
            enddo
         enddo
c     endif

      return
      end
!-----------------------------------------------------------------------
