!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParticle(imethod,ndim,npart,y)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer  imethod
      integer  ndim
      integer  npart
      real     y(*)

      call ppiclf_solve_InitDefault

      ppiclf_rparam(2)  = imethod
      ppiclf_rparam(12) = ndim
      ppiclf_npart      = npart

      call ppiclf_solve_InitParticleTag
      call ppiclf_solve_SetParticleTag

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
      call ppiclf_solve_SetupInterp

c     ! two-way coupling init
c     call ppiclf_project

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitDefault
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      ! set defaults
      ppiclf_rparam(2) = 0     ! time integ
c     ppiclf_rparam(3) =       ! poly order mesh... rm
c     ppiclf_rparam(4) =       ! 1=tracers, 0=projection... rm see param5

      ppiclf_rparam(5)  = -1   ! filt default for no projection

      ppiclf_rparam(8)  = 1    ! periodic in x (== 0) 
      ppiclf_rparam(9)  = 1    ! periodic in y (==0)
      ppiclf_rparam(10) = 1    ! periodic in z (==0)

      ppiclf_rparam(12) = 3    ! ndim


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

      ppiclf_rparam(5) = filt
      ppiclf_rparam(6) = real(ngrid)
      ppiclf_rparam(7) = alpha 

      rsig             = ppiclf_rparam(5)/(2.*sqrt(2.*log(2.)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_rparam(7)))

      call ppiclf_comm_CreateBin
      call ppiclf_comm_MapOverlapMesh

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
      subroutine ppiclf_solve_IntegrateParticle(time_,y,ydot)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      if (int(ppiclf_rparam(2)) .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(time_,y,ydot)

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
      call ppiclf_solve_SetRK3Coeff

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_user_SetYdot(time_,y,ydot)

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
      subroutine ppiclf_solve_SetupInterp
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

         call ppiclf_solve_RemoveParticle
         call ppiclf_comm_CreateBin
      if (ppiclf_rparam(5) .gt. 0) 
     >   call ppiclf_comm_MapOverlapMesh
         call ppiclf_comm_FindParticle
         call ppiclf_comm_MoveParticle

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(jp,fld)
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
      subroutine ppiclf_solve_ParallelProjection
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      if (int(ppiclf_rparam(4)) .ne. 1) then
         call ppiclf_comm_CreateGhost
         call ppiclf_comm_MoveGhost
         call ppiclf_solve_ProjectParticleGrid
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff
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
      subroutine ppiclf_solve_RemoveParticle
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
                    
         iproj(1,ip)  = ppiclf_iprop(8,ip)
         iproj(2,ip)  = ppiclf_iprop(9,ip)
         iproj(3,ip)  = ppiclf_iprop(10,ip)
         iproj(4,ip)  = ppiclf_iprop(11,ip)
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
                    
         iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
         iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
         iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
         iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp
c     ndum = ppiclf_npart_gp
c     ndum = ppiclf_npart

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
            if(ppiclf_rparam(12) .gt. 2) rdist2 = rdist2 +
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
      call fgslib_crystal_tuple_transfer(i_cr_hndl,neltbc,PPICLF_LEE
     >            , ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)

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
      subroutine ppiclf_solve_ProjectParticleBin
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
