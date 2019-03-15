!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParticle(imethod,ndim,npart,y)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer  imethod
      integer  ndim
      integer  npart
      real     y(*)

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitParticle$',0.
     >   ,ppiclf_nid)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter must be before InitParticle$',0.,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.,0)

      call ppiclf_prints('*Begin InitParticle$')

         call ppiclf_prints('   *Begin InitParam$')
            call ppiclf_solve_InitParam(imethod,ndim,npart)
         call ppiclf_prints('    End InitParam$')
         
            call ppiclf_copy(ppiclf_y,y,ppiclf_npart)
         
         call ppiclf_prints('   *Begin ParticleTag$')
            call ppiclf_solve_InitParticleTag
            call ppiclf_solve_SetParticleTag
         call ppiclf_prints('    End ParticleTag$')
         
c        call ppiclf_prints('   *Begin RemoveParticle$')
c           call ppiclf_solve_RemoveParticle
c        call ppiclf_prints('    End RemoveParticle$')
         
         call ppiclf_prints('   *Begin CreateBin$')
            call ppiclf_comm_CreateBin
         call ppiclf_prints('    End CreateBin$')
         
         call ppiclf_prints('   *Begin FindParticle$')
            call ppiclf_comm_FindParticle
         call ppiclf_prints('    End FindParticle$')
         
         call ppiclf_prints('   *Begin MoveParticle$')
            call ppiclf_comm_MoveParticle
         call ppiclf_prints('    End MoveParticle$')

            call ppiclf_solve_OutputDiagGen

      call ppiclf_prints(' End InitParticle$')

      PPICLF_LINIT = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_OutputDiagAll
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_solve_OutputDiagGen
      call ppiclf_solve_OutputDiagGhost
      if (ppiclf_lfilt) call ppiclf_solve_OutputDiagSubBin
      if (ppiclf_overlap) call ppiclf_solve_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_OutputDiagGen
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = PPICLF_NDXGP*PPICLF_NDYGP*PPICLF_NDZGP

      call ppiclf_printsi('  -Cycle                  :$',ppiclf_cycle)
      call ppiclf_printsi('  -Output Freq.           :$',ppiclf_iostep)
      call ppiclf_printsr('  -Time                   :$',ppiclf_time)
      call ppiclf_printsr('  -dt                     :$',ppiclf_dt)
      call ppiclf_printsi('  -Global particles       :$',npart_tot)
      call ppiclf_printsi('  -Local particles (Max)  :$',npart_max)
      call ppiclf_printsi('  -Local particles (Min)  :$',npart_min)
      call ppiclf_printsi('  -Local particles (Ideal):$',npart_ide)
      call ppiclf_printsi('  -Total ranks            :$',ppiclf_np)
      call ppiclf_printsi('  -Problem dimensions     :$',ppiclf_ndim)
      call ppiclf_printsi('  -Integration method     :$',ppiclf_imethod)
      call ppiclf_printsi('  -Number of bins total   :$',nbin_total)
      call ppiclf_printsi('  -Number of bins (x)     :$',PPICLF_NDXGP)
      call ppiclf_printsi('  -Number of bins (y)     :$',PPICLF_NDYGP)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$',PPICLF_NDZGP)
      call ppiclf_printsr('  -Bin xl coordinate      :$',ppiclf_binb(1))
      call ppiclf_printsr('  -Bin xr coordinate      :$',ppiclf_binb(2))
      call ppiclf_printsr('  -Bin yl coordinate      :$',ppiclf_binb(3))
      call ppiclf_printsr('  -Bin yr coordinate      :$',ppiclf_binb(4))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zl coordinate      :$',ppiclf_binb(5))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zr coordinate      :$',ppiclf_binb(6))

      call ppiclf_prints('  End General Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_OutputDiagGrid
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_prints(' *Begin Grid Info$')

         nel_max_orig   = ppiclf_iglmax(ppiclf_nee,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nee,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nee,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_neltb,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_neltb,1)
         nel_total_map = ppiclf_iglsum(ppiclf_neltb,1)

      call ppiclf_printsi('  -Orig. Global cells     :$',nel_total_orig)
      call ppiclf_printsi('  -Orig. Local cells (Max):$',nel_max_orig)
      call ppiclf_printsi('  -Orig. Local cells (Min):$',nel_min_orig)
      call ppiclf_printsi('  -Map Global cells       :$',nel_total_map)
      call ppiclf_printsi('  -Map Local cells (Max)  :$',nel_max_map)
      call ppiclf_printsi('  -Map Local cells (Min)  :$',nel_min_map)

      call ppiclf_prints('  End Grid Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_OutputDiagGhost
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_prints(' *Begin Ghost Info$')

         npart_max = ppiclf_iglmax(ppiclf_npart_gp,1)
         npart_min = ppiclf_iglmin(ppiclf_npart_gp,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart_gp,1)

      call ppiclf_printsi('  -Global ghosts          :$',npart_tot)
      call ppiclf_printsi('  -Local ghosts (Max)     :$',npart_max)
      call ppiclf_printsi('  -Local ghosts (Min)     :$',npart_min)

      call ppiclf_prints('  End Ghost Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_OutputDiagSubBin
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      call ppiclf_prints(' *Begin SubBin Info$')

      nbin_total = ppiclf_bx*ppiclf_by*ppiclf_bz

      call ppiclf_printsi('  -Number of local bin    :$',nbin_total)
      call ppiclf_printsi('  -Number of local bin (x):$',PPICLF_BX)
      call ppiclf_printsi('  -Number of local bin (y):$',PPICLF_BY)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of local bin (z):$',PPICLF_BZ)
      call ppiclf_printsr('  -Bin width (x)          :$',ppiclf_rdx)
      call ppiclf_printsr('  -Bin width (y)          :$',ppiclf_rdy)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin width (z)          :$',ppiclf_rdz)
      call ppiclf_printsr('  -Filter width           :$',ppiclf_filter)
      call ppiclf_printsr('  -Filter cut-off         :$'
     >                                                 ,ppiclf_d2chk(2))
      call ppiclf_printsi('  -SubBins per filter res.:$',ppiclf_ngrids)

      call ppiclf_prints('  End SubBin Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParam(imethod,ndim,npart)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      integer  imethod
      integer  ndim
      integer  npart

      if (imethod .le. 0 .or. imethod .ge. 2)
     >   call ppiclf_exittr('Invalid integration method$',0.0,imethod)
      if (ndim .le. 1 .or. ndim .ge. 4)
     >   call ppiclf_exittr('Invalid problem dimension$',0.0,ndim)
      if (npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   call ppiclf_exittr('Invalid number of particles$',0.0,npart)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_npart        = npart

      ppiclf_filter       = -1   ! filt default for no projection

      ppiclf_iperiodic(1) = 1    ! periodic in x (== 0) 
      ppiclf_iperiodic(2) = 1    ! periodic in y (==0)
      ppiclf_iperiodic(3) = 1    ! periodic in z (==0)

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0
      ppiclf_time   = 0.0

      ppiclf_restart = .false.
      ppiclf_overlap = .false.
      ppiclf_linit   = .false.
      ppiclf_lfilt   = .false.

      ppiclf_xdrange(1,1) = -1E8
      ppiclf_xdrange(2,1) =  1E8
      ppiclf_xdrange(1,2) = -1E8
      ppiclf_xdrange(2,2) =  1E8
      ppiclf_xdrange(1,3) = -1E8
      ppiclf_xdrange(2,3) =  1E8

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,ngrid)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real    filt
      real    alpha
      integer ngrid

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.,0)

      ppiclf_filter = filt
      ppiclf_ngrids = ngrid
      ppiclf_alpha  = alpha 

      rsig             = ppiclf_filter/(2.*sqrt(2.*log(2.)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_alpha))

      call ppiclf_prints('   *Redo CreateBin$')
         call ppiclf_comm_CreateBin
      call ppiclf_prints('    End CreateBin$')

      call ppiclf_prints('   *Begin CreateSubBin$')
         call ppiclf_comm_CreateSubBin
      call ppiclf_prints('    End CreateSubBin$')

      call ppiclf_prints('   *Redo FindParticle$')
         call ppiclf_comm_FindParticle
      call ppiclf_prints('    End FindParticle$')

      call ppiclf_prints('   *Redo MoveParticle$')
         call ppiclf_comm_MoveParticle
      call ppiclf_prints('    End MoveParticle$')

      call ppiclf_solve_OutputDiagSubBin

      PPICLF_LFILT = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetParticleTag
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      do i=1,ppiclf_npart
         if (ppiclf_iprop(5,i) .eq. -1) ppiclf_iprop(5,i) = ppiclf_nid 
         if (ppiclf_iprop(6,i) .eq. -1) ppiclf_iprop(6,i) = ppiclf_cycle
         if (ppiclf_iprop(7,i) .eq. -1) ppiclf_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParticleTag
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
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time
     >                                         ,y,ydot)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer istep
      integer iostep
      real    dt
      real    time
      real    y(*)
      real    ydot(*)

      integer icalld
      save    icalld
      data    icalld /0/

      ! write initial condition (with projection)
      if (icalld .eq. 0) then
         icalld = icalld + 1

         call ppiclf_solve_OutputDiagAll

         call ppiclf_io_WriteParticleVTU('',0)

         if (ppiclf_lfilt)
     >      call ppiclf_io_WriteSubBinVTU('',0)
      endif


      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      if (ppiclf_imethod .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(time,y,ydot)


      ! output files
      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0) then

         call ppiclf_solve_OutputDiagAll

         call ppiclf_io_WriteParticleVTU('',0)

         if (ppiclf_lfilt)
     >      call ppiclf_io_WriteSubBinVTU('',0)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3(time_,y,ydot)
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
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_user_SetYdot(time_,y,ydot)

         ndum = PPICLF_NPART*PPICLF_LRS
         ! rk3 integrate
         do i=1,ndum
            y(i) =  ppiclf_rk3coef(1,istage)*ppiclf_y1 (i)
     >            + ppiclf_rk3coef(2,istage)*y      (i)
     >            + ppiclf_rk3coef(3,istage)*ydot   (i)
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSolve
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

c        call ppiclf_solve_RemoveParticle
         call ppiclf_comm_CreateBin
         if (ppiclf_lfilt) call ppiclf_comm_CreateSubBin
         if (ppiclf_overlap) call ppiclf_comm_MapOverlapMesh
         call ppiclf_comm_FindParticle
         call ppiclf_comm_MoveParticle

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(jp,fld)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real fld(*)
      logical partl

      ! should group all fields together if more than one ...

      ! also throw error if overlap is not set

c     real xm1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE), 
c    >     ym1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE), 
c    >     zm1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
c     common /ppiclf_tmp_grid/ xm1, ym1, zm1

c     n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
c     do ie=1,ppiclf_neltb
c        call ppiclf_copy(xm1(1,1,1,ie),ppiclf_xm1b(1,1,1,1,ie),n)
c        call ppiclf_copy(ym1(1,1,1,ie),ppiclf_xm1b(1,1,1,2,ie),n)
c        call ppiclf_copy(zm1(1,1,1,ie),ppiclf_xm1b(1,1,1,3,ie),n)
c     enddo

c     tol     = 5e-13
c     bb_t    = 0.01
c     npt_max = 128

c     call fgslib_findpts_setup(ppiclf_fp_hndl
c    >                         ,ppiclf_comm_nid
c    >                         ,1 ! only 1 rank on this comm
c    >                         ,ppiclf_ndim
c    >                         ,xm1
c    >                         ,ym1
c    >                         ,zm1
c    >                         ,PPICLF_LEX
c    >                         ,PPICLF_LEY
c    >                         ,PPICLF_LEZ
c    >                         ,ppiclf_neltb
c    >                         ,2*PPICLF_LEX
c    >                         ,2*PPICLF_LEY
c    >                         ,2*PPICLF_LEZ
c    >                         ,bb_t
c    >                         ,ppiclf_neltb+2
c    >                         ,ppiclf_neltb+2
c    >                         ,npt_max
c    >                         ,tol)


c     neltbc = ppiclf_neltbb
c     do ie=1,neltbc
c        call ppiclf_icopy(ppiclf_er_mapc(1,ie),ppiclf_er_maps(1,ie)
c    >             ,PPICLF_LRMAX)
c     enddo

c     nl   = 0
c     nii  = PPICLF_LRMAX
c     njj  = 6
c     nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
c     nrr  = nxyz*1
c     nkey = 3
c     call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc
c    >      ,PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,fld,nrr,njj)
c     call fgslib_crystal_tuple_sort    (ppiclf_cr_hndl,neltbc
c    >       ,ppiclf_er_mapc,nii,partl,nl,fld,nrr,nkey,1)

c     ix = 1
c     iy = 2
c     iz = 3

c     call fgslib_findpts(PPICLF_FP_HNDL           !   call fgslib_findpts( ihndl,
c    $        , ppiclf_iprop (1 ,1),PPICLF_LIP        !   $             rcode,1,
c    $        , ppiclf_iprop (3 ,1),PPICLF_LIP        !   &             proc,1,
c    $        , ppiclf_iprop (2 ,1),PPICLF_LIP        !   &             elid,1,
c    $        , ppiclf_rprop2(1 ,1),PPICLF_LRP2       !   &             rst,ndim,
c    $        , ppiclf_rprop2(4 ,1),PPICLF_LRP2       !   &             dist,1,
c    $        , ppiclf_y     (ix,1),PPICLF_LRS        !   &             pts(    1),1,
c    $        , ppiclf_y     (iy,1),PPICLF_LRS        !   &             pts(  n+1),1,
c    $        , ppiclf_y     (iz,1),PPICLF_LRS ,PPICLF_NPART) !   &             pts(2*n+1),1,n)

c     call fgslib_findpts_eval_local( PPICLF_FP_HNDL
c    >                               ,ppiclf_rprop (jp,1)
c    >                               ,PPICLF_LRP
c    >                               ,ppiclf_iprop (2,1)
c    >                               ,PPICLF_LIP
c    >                               ,ppiclf_rprop2(1,1)
c    >                               ,PPICLF_LRP2
c    >                               ,PPICLF_NPART
c    >                               ,fld)

c     call fgslib_findpts_free(PPICLF_FP_HNDL)

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ParallelProjection
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      if (ppiclf_lfilt) then
         call ppiclf_comm_CreateGhost
         call ppiclf_comm_MoveGhost
         if (ppiclf_overlap) call ppiclf_solve_ProjectParticleGrid
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff(dt)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real dt

      ppiclf_rk3coef(1,1) = 0.0
      ppiclf_rk3coef(2,1) = 1.0 
      ppiclf_rk3coef(3,1) = dt
      ppiclf_rk3coef(1,2) = 3.0/4.0
      ppiclf_rk3coef(2,2) = 1.0/4.0 
      ppiclf_rk3coef(3,2) = dt/4.0 
      ppiclf_rk3coef(1,3) = 1.0/3.0
      ppiclf_rk3coef(2,3) = 2.0/3.0 
      ppiclf_rk3coef(3,3) = dt*2.0/3.0 


      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_RemoveParticle
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >                                   iperiodicz,ndim

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
      ndim       = ppiclf_ndim

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
               call ppiclf_copy
     >              (ppiclf_y     (1,ic),ppiclf_y(1,i)     ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_y1    (isr) ,ppiclf_y1(isl)    ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydot  (1,ic),ppiclf_ydot(1,i)  ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydotc (1,ic),ppiclf_ydotc(1,i) ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_rprop (1,ic),ppiclf_rprop(1,i) ,PPICLF_LRP)
               call ppiclf_copy
     >              (ppiclf_rprop2(1,ic),ppiclf_rprop2(1,i),PPICLF_LRP2)
               call ppiclf_icopy
     >              (ppiclf_iprop(1,ic) ,ppiclf_iprop(1,i) ,PPICLF_LIP)
            endif
         endif
      enddo
      ppiclf_npart = ic

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleGrid
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real    multfci

      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      integer ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp

      logical partl

      real pi

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      nxyzdum = nxyz*PPICLF_LRP_PRO*PPICLF_LEE
      do i=1,nxyzdum
         ppiclf_pro_fldb(i,1,1,1,1) = 0.0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 3

      ! real particles
      do ip=1,ppiclf_npart
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_ndim .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip) = rbexpi
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)


         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip) = ppiclf_cp_map(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip)  = ppiclf_iprop(8,ip)
         iproj(2,ip)  = ppiclf_iprop(9,ip)
         iproj(3,ip)  = ppiclf_iprop(10,ip)
         iproj(4,ip)  = ppiclf_iprop(11,ip)
      enddo

      ! ghost particles
      do ip=1,ppiclf_npart_gp
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_ndim .gt. 2) multfci = multfci**(1.5d+0)
         rbexpi   = 1./(-2.*rsig**2)

         rproj(1 ,ip+ppiclf_npart) = rbexpi
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         do j=5,PPICLF_LRP_GP+1
            rproj(j,ip+ppiclf_npart) = ppiclf_rprop_gp(j-1,ip)*multfci
         enddo
                    
         iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
         iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
         iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
         iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp

      do ip=1,ndum
         iip      = iproj(1,ip)
         jjp      = iproj(2,ip)
         kkp      = iproj(3,ip)
         ndumdum  = iproj(4,ip)

         ilow  = iip-1
         ihigh = iip+1
         jlow  = jjp-1
         jhigh = jjp+1
         klow  = kkp-1
         khigh = kkp+1

         do ie=1,ppiclf_neltb

               if (ppiclf_el_map(1,ie) .gt. ndumdum) exit
               if (ppiclf_el_map(2,ie) .lt. ndumdum) cycle 
         
               if (ppiclf_el_map(3,ie) .gt. ihigh) cycle
               if (ppiclf_el_map(4,ie) .lt. ilow)  cycle
               if (ppiclf_el_map(5,ie) .gt. jhigh) cycle
               if (ppiclf_el_map(6,ie) .lt. jlow)  cycle
               if (ppiclf_el_map(7,ie) .gt. khigh) cycle
               if (ppiclf_el_map(8,ie) .lt. klow)  cycle

         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
            if (ppiclf_modgp(i,j,k,ie,4).ne.ndumdum) cycle

            rdist2  = (ppiclf_xm1b(i,j,k,1,ie) - rproj(2,ip))**2 +
     >                (ppiclf_xm1b(i,j,k,2,ie) - rproj(3,ip))**2
            if(ppiclf_ndim .gt. 2) rdist2 = rdist2 +
     >                (ppiclf_xm1b(i,j,k,3,ie) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = exp(rdist2*rproj(1,ip))
            
            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                         ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                       + rproj(j1,ip)*rexp
            enddo
         enddo
         enddo
         enddo
         enddo
      enddo

      ! now send xm1b to the processors in nek that hold xm1

      neltbc = ppiclf_neltb
      ndum = PPICLF_LRMAX*neltbc
      call ppiclf_icopy(ppiclf_er_mapc,ppiclf_er_map,ndum)
      do ie=1,neltbc
         ppiclf_er_mapc(5,ie) = ppiclf_er_mapc(2,ie)
         ppiclf_er_mapc(6,ie) = ppiclf_er_mapc(2,ie)
      enddo
      nl = 0
      nii = PPICLF_LRMAX
      njj = 6
      nrr = nxyz*PPICLF_LRP_PRO
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc,
     >   PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)

      ! add the fields from the bins to ptw array
      nlxyzep = nxyz*PPICLF_LEE*PPICLF_LRP_PRO
      do i=1,nlxyzep
         PPICLF_PRO_FLD(i,1,1,1,1) = 0.0
      enddo


      do ie=1,neltbc
         iee = ppiclf_er_mapc(1,ie)
         do ip=1,PPICLF_LRP_PRO
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
           PPICLF_PRO_FLD(i,j,k,iee,ip) = PPICLF_PRO_FLD(i,j,k,iee,ip) +
     >                                    PPICLF_PRO_FLDB(i,j,k,ip,ie)
         enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleSubBin
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real    multfci
      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      integer ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp

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

      ! real particles
      do ip=1,ppiclf_npart
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_ndim .gt. 2) multfci = multfci**(1.5d+0)
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
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (ppiclf_ndim .gt. 2) multfci = multfci**(1.5d+0)
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

      idum = floor(ppiclf_filter/2.0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1
      jdum = floor(ppiclf_filter/2.0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1
      kdum = floor(ppiclf_filter/2.0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1

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

         do k=kl,kr
         do j=jl,jr
         do i=il,ir
    
            rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
            if(ppiclf_ndim .gt. 2) rdist2 = rdist2 +
     >                (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = exp(rdist2*rproj(1,ip))
            
            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_grid_fld(i,j,k,jj) = 
     >                         ppiclf_grid_fld(i,j,k,jj) 
     >                       + rproj(j1,ip)*rexp
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
