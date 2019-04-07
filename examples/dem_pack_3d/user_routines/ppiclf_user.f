!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(time,y,ydot)
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Input:
!
      real*8 time
      real*8 y(*)
!
! Output:
!
      real*8 ydot(*)
!
! Internal:
!
      real*8 rmass,fbx,fby,fbz,fcx,fcy,fcz
      integer*4 i,j
!
! evaluate ydot
      do i=1,ppiclf_npart
         ! striding solution y vector
         j = PPICLF_LRS*(i-1)

         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Gravity
         fbx  = 0.0
         fby  = -9.8*rmass
         fbz  = 0.0

         call ppiclf_solve_NearestNeighbor(i)

         fcx  = ppiclf_ydotc(PPICLF_JVX,i)
         fcy  = ppiclf_ydotc(PPICLF_JVY,i)
         fcz  = ppiclf_ydotc(PPICLF_JVZ,i)

         ! set ydot for all PPICLF_SLN number of equations
         ydot(PPICLF_JX +j) = y(PPICLF_JVX +j)
         ydot(PPICLF_JY +j) = y(PPICLF_JVY +j)
         ydot(PPICLF_JZ +j) = y(PPICLF_JVZ +j)
         ydot(PPICLF_JVX+j) = (fbx+fcx)/rmass
         ydot(PPICLF_JVY+j) = (fby+fcy)/rmass
         ydot(PPICLF_JVZ+j) = (fbz+fcz)/rmass
      enddo 
c evaluate ydot

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
! Input:
!
      real*8 y(*)
      real*8 ydot(*)
      real*8 ydotc(*)
      real*8 rprop(*)
!
! Output:
!
      real*8 map(*)
!

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Input:
!
      integer*4 i
      integer*4 j
      real*8 yi(*)     ! PPICLF_LRS
      real*8 rpropi(*) ! PPICLF_LRP
      real*8 yj(*)     ! PPICLF_LRS
      real*8 rpropj(*) ! PPICLF_LRP
!
! Internal:
!
      ! For user implemented collision model
      real*8 ksp,erest
      common /external_user_collsion/ ksp,erest
      ! For user implemented collision model
      real*8 rpi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra
!
      rpi2  =  9.869604401089358

      ! other particles
      if (j .ne. 0) then
         rthresh  = 0.5*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2 + rzdiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
         rmult = 1./sqrt(1./rm1+1./rm2)
         eta   = 2.*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x +
     >           (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y +
     >           (yj(PPICLF_JVZ)-yi(PPICLF_JVZ))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         PPICLF_YDOTC(PPICLF_JVZ,i) = PPICLF_YDOTC(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z

      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp

         ! give a bit larger collision threshold for walls
         rextra   = 0.5
         rthresh  = (0.5+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2 + rzdiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -1.0*(yi(PPICLF_JVX)*rn_12x +
     >                    yi(PPICLF_JVY)*rn_12y +
     >                    yi(PPICLF_JVZ)*rn_12z)

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         PPICLF_YDOTC(PPICLF_JVZ,i) = PPICLF_YDOTC(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z

      endif

      return
      end
!-----------------------------------------------------------------------
