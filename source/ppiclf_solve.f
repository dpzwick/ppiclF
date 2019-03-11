!-----------------------------------------------------------------------
      subroutine ppiclf_init(rparam,y,npart,time_)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real     y(*)
      real     rparam(*)

      ppiclf_d2chk(2) = 0.0
      ppiclf_npart = npart

      call ppiclf_comm_setup

      call ppiclf_rparam_set(rparam)
      call ppiclf_tag_init
      call ppiclf_tag_set
      call ppiclf_init_filter

c     ! get domain bounds
      ppiclf_xdrange(1,1) = -1E8
      ppiclf_xdrange(2,1) =  1E8
      ppiclf_xdrange(1,2) = -1E8
      ppiclf_xdrange(2,2) =  1E8
      ppiclf_xdrange(1,3) = -1E8
      ppiclf_xdrange(2,3) =  1E8
c     call domain_size( ppiclf_xdrange(1,1),ppiclf_xdrange(2,1)
c    >                 ,ppiclf_xdrange(1,2),ppiclf_xdrange(2,2)
c    >                 ,ppiclf_xdrange(1,3),ppiclf_xdrange(2,3))

      ppiclf_nbmax = PPICLF_BMAX ! set only one bin for now...

c     ! send particles to correct rank
      call ppiclf_interpolate_setup

c     ! two-way coupling init
      call ppiclf_project

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_rparam_set(rparam)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real rparam(*)

c     call rzero(ppiclf_rparam, ppiclf_nparam)
      do i=1,ppiclf_nparam
         ppiclf_rparam(i) = 0.0
      enddo

      ! set defaults
      if (rparam(1) .eq. 0) then
         ppiclf_rparam(1)  = 0        ! use custom values
         ppiclf_rparam(2)  = 1        ! time integration method
         ppiclf_rparam(3)  = PPICLF_LX1  ! polynomial order of mesh
         ppiclf_rparam(4)  = 1        ! use 1 for tracers only
         ppiclf_rparam(5)  = 0        ! filter size in real units

         ppiclf_rparam(7)  = 0        ! percent decay of Gaussian filter
         ppiclf_rparam(8)  = 1        ! periodic in x
         ppiclf_rparam(9)  = 1        ! periodic in y
         ppiclf_rparam(10) = 1        ! periodic in z
         ppiclf_rparam(11) = -1       ! time step
         ppiclf_rparam(12) = 3        ! problem dimension

      ! custom values
      else
         do i=1,ppiclf_nparam
            ppiclf_rparam(i) = rparam(i)
         enddo

         tol = 1.e-12
         if (wdsize.eq.4) tol = 1.e-6

         if (abs(ppiclf_rparam(2)) .lt. tol  ) 
     >      ppiclf_rparam(2) = 1
         if (abs(ppiclf_rparam(3)) .lt. tol  ) 
     >      ppiclf_rparam(3) = PPICLF_LX1

      endif

      ppiclf_dt = ppiclf_rparam(11)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_init_filter
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      rsig  = 0.0
      alph  = ppiclf_rparam(7)

      rsig  = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
      ppiclf_d2chk(2) = rsig*sqrt(-2*log(alph))

c     if (int(ppiclf_rparam(4)) .eq. 1) ppiclf_d2chk(2) = 0

c     rdx_max = 0.0

c     lx1m1 = max(1,PPICLF_LX1-1)
c     ly1m1 = max(1,PPICLF_LY1-1)
c     lz1m1 = max(1,PPICLF_LZ1-1)

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

c     ppiclf_rdx_max = rdx_max

c     rsig_dx_min_set  = 0.25
c     rfilt_dp_min_set = 1.0
c     rsig_dx_min      = 1E12

c     tol = 1.e-12
c     if (wdsize.eq.4) tol = 1.e-6

c     rfilt = ppiclf_rparam(6)
c     ppiclf_d2chk(2) = -1
c     if (abs(alph) .lt. tol) then
c        alph = 1E-3
c        ppiclf_rparam(7) = alph
c     endif

c     do i=1,ppiclf_npart

c        if (filt .lt. tol) then
c           rsig_test_grid = rsig_dx_min_set*rdx_max + tol
c           rsig_test_diam = rfilt_dp_min_set*ppiclf_rprop(jdp,i)
c    >                       /(2.*sqrt(2.*log(2.)))
c           rsig_test = max(rsig_test_grid,rsig_test_diam)
c           filt      = rsig_test/ppiclf_rprop(jdp,i)*2.*sqrt(2.*log(2.))
c           rfilt = filt
c        else
c           rsig_test = filt*ppiclf_rprop(jdp,i)/(2.*sqrt(2.*log(2.)))
c        endif

c        rsig_dx = rsig_test/rdx_max
c        if (rsig_dx .lt. rsig_dx_min) rsig_dx_min = rsig_dx

c        if (rsig_test .gt. rsig) rsig = rsig_test
c     enddo
c     rsig = glmax(rsig,1)
c     rsig_dx_min = glmin(rsig_dx_min,1)
c        
c     ppiclf_d2chk(2) = rsig*sqrt(-2*log(alph))

c     if (rsig_dx_min .lt. rsig_dx_min_set) then
c        if (nid .eq. 0) then
c           filt_new = filt/(rsig_dx_min/rsig_dx_min_set)
c           write(6,100) filt, filt_new
c100        format('Reset Gaussian filter width from', E14.7
c    >             ' to', E14.7, ' or larger')
c           call exitt
c        endif
c     endif

c     ppiclf_rparam(6) = glmax(rfilt,1)
c     ppiclf_d2chk(2)  = glmax(ppiclf_d2chk(2),1)

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_tag_set
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      do i=1,ppiclf_npart
         if (ppiclf_iprop(5,i) .eq. -1) ppiclf_iprop(5,i) = 0 ! nid
         if (ppiclf_iprop(6,i) .eq. -1) ppiclf_iprop(6,i) = 0 ! istep
         if (ppiclf_iprop(7,i) .eq. -1) ppiclf_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_tag_init
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      if (.not. PPICLF_RESTART) then
         do i=1,PPICLF_LPART
            ppiclf_iprop(5,i) = -1
            ppiclf_iprop(6,i) = -1
            ppiclf_iprop(7,i) = -1
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve(time_,y,ydot)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      if (int(ppiclf_rparam(2)) .eq. 1) 
     >   call ppiclf_rk3_driver(time_,y,ydot)

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_rk3_driver(time_,y,ydot)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      ndum = PPICLF_NPART*PPICLF_LRS

      ! save stage 1 solution
      do i=1,ndum
         ppiclf_y1(i) = y(i)
      enddo

      ! get rk3 coeffs
      call ppiclf_rk3_coeff

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_fun(time_,y,ydot)

         ndum = PPICLF_NPART*PPICLF_LRS
         ! rk3 integrate
         do i=1,ndum
            y(i) =  rk3coef(1,istage)*ppiclf_y1 (i)
     >            + rk3coef(2,istage)*y      (i)
     >            + rk3coef(3,istage)*ydot   (i)
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_interpolate_setup
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_move_outlier
      call ppiclf_comm_bin_setup
      call ppiclf_comm_findpts
      call ppiclf_comm_crystal

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_interpolate_fld(jp,fld)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      common /intp_h/ ih_intp(2,1)

      real fld(*)

      ih_intp1 = ih_intp(1,i_fp_hndl)

c     call fgslib_findpts_eval_local( ih_intp1
c    >                               ,ppiclf_rprop (jp,1)
c    >                               ,PPICLF_LRP
c    >                               ,ppiclf_iprop (2,1)
c    >                               ,PPICLF_LIP
c    >                               ,ppiclf_rprop2(1,1)
c    >                               ,PPICLF_LRP2
c    >                               ,PPICLF_NPART
c    >                               ,fld)

      ! this one
c        call fgslib_findpts_eval( ih_intp1
c    >                           ,ppiclf_rprop (jp,1)
c    >                           ,PPICLF_LRP
c    >                           ,ppiclf_iprop (1,1)
c    >                           ,PPICLF_LIP
c    >                           ,ppiclf_iprop (3,1)
c    >                           ,PPICLF_LIP
c    >                           ,ppiclf_iprop (2,1)
c    >                           ,PPICLF_LIP
c    >                           ,ppiclf_rprop2(1,1)
c    >                           ,PPICLF_LRP2
c    >                           ,PPICLF_NPART
c    >                           ,fld)

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_project
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      if (int(ppiclf_rparam(4)) .ne. 1) then
         call ppiclf_comm_ghost_create
         call ppiclf_comm_ghost_send
         call ppiclf_solve_project
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_rk3_coeff
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      rk3coef(1,1) = 0.0
      rk3coef(2,1) = 1.0 
      rk3coef(3,1) = ppiclf_dt
      rk3coef(1,2) = 3.0/4.0
      rk3coef(2,2) = 1.0/4.0 
      rk3coef(3,2) = ppiclf_dt/4.0 
      rk3coef(1,3) = 1.0/3.0
      rk3coef(2,3) = 2.0/3.0 
      rk3coef(3,3) = ppiclf_dt*2.0/3.0 

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_move_outlier
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >                                   iperiodicz

      iperiodicx = int(ppiclf_rparam(8))
      iperiodicy = int(ppiclf_rparam(9))
      iperiodicz = int(ppiclf_rparam(10))
      ndim       = int(ppiclf_rparam(12))

      jj(1) = 1
      jj(2) = 2
      jj(3) = 3

      do i=1,ppiclf_npart
         isl = (i -1) * PPICLF_LRS + 1
         in_part(i) = 0
         do j=0,ndim-1
            jchk = jj(j+1)
            if (ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
                   ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
                  goto 1512
                endif
            endif
            if (ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) +
     &                     abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
                   ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
     &                     abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (ppiclf_iprop(1,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1512 continue
         enddo
      enddo

      ic = 0
      do i=1,ppiclf_npart
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               isl = (i -1) * PPICLF_LRS + 1
               isr = (ic-1) * PPICLF_LRS + 1
               call copy
     >              (ppiclf_y     (1,ic),ppiclf_y(1,i)     ,PPICLF_LRS)
               call copy
     >              (ppiclf_y1    (isr) ,ppiclf_y1(isl)    ,PPICLF_LRS)
               call copy
     >              (ppiclf_ydot  (1,ic),ppiclf_ydot(1,i)  ,PPICLF_LRS)
               call copy
     >              (ppiclf_ydotc (1,ic),ppiclf_ydotc(1,i) ,PPICLF_LRS)
               call copy
     >              (ppiclf_rprop (1,ic),ppiclf_rprop(1,i) ,PPICLF_LRP)
               call copy
     >              (ppiclf_rprop2(1,ic),ppiclf_rprop2(1,i),PPICLF_LRP2)
               call icopy
     >              (ppiclf_iprop(1,ic) ,ppiclf_iprop(1,i) ,PPICLF_LIP)
            endif
         endif
      enddo
      ppiclf_npart = ic

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_project
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real    multfci
      integer e

      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      logical partl

      real pi

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 3

c     ppiclf_npart_gp = 0

      ! real particles
      do ip=1,ppiclf_npart
         rsig    = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip) = rbexpi
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)


         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip) = ppiclf_cp_map(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - ppiclf_biny(1,1))/ppiclf_rdy)
         iproj(3,ip) = 
     >       floor( (rproj(4,ip) - ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ! ghost particles
      do ip=1,ppiclf_npart_gp
         rsig    = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip+ppiclf_npart) = rbexpi
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip+ppiclf_npart) = ppiclf_rprop_gp(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip+ppiclf_npart) = 
     >     floor((rproj(2,ip+ppiclf_npart)-ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip+ppiclf_npart) = 
     >     floor((rproj(3,ip+ppiclf_npart)-ppiclf_biny(1,1))/ppiclf_rdy)
         iproj(3,ip+ppiclf_npart) = 
     >     floor((rproj(4,ip+ppiclf_npart)-ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp
c     ndum = ppiclf_npart_gp
c     ndum = ppiclf_npart

      idum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1
      jdum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1
      kdum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1

      do ip=1,ndum
      ! project here
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_project_bins
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real    multfci
      integer e

      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      logical partl

      real pi

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 3

c     ppiclf_npart_gp = 0

      ! real particles
      do ip=1,ppiclf_npart
         rsig    = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip) = rbexpi
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)


         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip) = ppiclf_cp_map(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - ppiclf_biny(1,1))/ppiclf_rdy)
         iproj(3,ip) = 
     >       floor( (rproj(4,ip) - ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ! ghost particles
      do ip=1,ppiclf_npart_gp
         rsig    = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_rparam(12) .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip+ppiclf_npart) = rbexpi
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip+ppiclf_npart) = ppiclf_rprop_gp(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip+ppiclf_npart) = 
     >     floor((rproj(2,ip+ppiclf_npart)-ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip+ppiclf_npart) = 
     >     floor((rproj(3,ip+ppiclf_npart)-ppiclf_biny(1,1))/ppiclf_rdy)
         iproj(3,ip+ppiclf_npart) = 
     >     floor((rproj(4,ip+ppiclf_npart)-ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp
c     ndum = ppiclf_npart_gp
c     ndum = ppiclf_npart

      idum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1
      jdum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1
      kdum = floor(ppiclf_rparam(5)/2.0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_rparam(7))/log(2.0)))+1

c     if (int(ppiclf_rparam(3)) .eq. 1) then
         do ip=1,ndum
            iip = iproj(1,ip)
            jjp = iproj(2,ip)
            kkp = iproj(3,ip)

            il  = max(1     ,iip-idum)
            ir  = min(ppiclf_bx,iip+idum)
            jl  = max(1     ,jjp-jdum)
            jr  = min(ppiclf_by,jjp+jdum)
            kl  = max(1     ,kkp-kdum)
            kr  = min(ppiclf_bz,kkp+kdum)

c           write(6,*) il,ir,jl,jr,kl,kr

c           do k=1,ppiclf_bz
c           do j=1,ppiclf_by
c           do i=1,ppiclf_bx
            do k=kl,kr
            do j=jl,jr
            do i=il,ir
    
c              if (ppiclf_grid_ii(i,j,k) .gt. ir) cycle
c              if (ppiclf_grid_ii(i,j,k) .lt. il) cycle
c              if (ppiclf_grid_jj(i,j,k) .gt. jr) cycle
c              if (ppiclf_grid_jj(i,j,k) .lt. jl) cycle
c              if (ppiclf_grid_kk(i,j,k) .gt. kr) cycle
c              if (ppiclf_grid_kk(i,j,k) .lt. kl) cycle


               rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                   (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
               if(ppiclf_rparam(12) .gt. 2) rdist2 = rdist2 +
     >                   (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

               if (rdist2 .gt. d2chk2_sq) cycle

c              write(6,*) 'Made itt', ppiclf_grid_x(i,j,k)
c    >                              , ppiclf_grid_y(i,j,k)
c    >                              , ppiclf_grid_z(i,j,k)

            
               rexp = exp(rdist2*rproj(1,ip))
               
               do jj=1,PPICLF_LRP_PRO
                  j1 = jj+4
                  ppiclf_grid_fld(i,j,k,jj) = 
     >                            ppiclf_grid_fld(i,j,k,jj) 
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
