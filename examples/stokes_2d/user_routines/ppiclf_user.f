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
      real*8 fqsx,fqsy,fbx,fby
      integer*4 i,j
!
! evaluate ydot
      do i=1,ppiclf_npart
         ! striding solution y vector
         j = PPICLF_LRS*(i-1)

         ! Stokes drag
         fqsx = -y(PPICLF_JVX+j)/ppiclf_rprop(PPICLF_R_JTAUP,i)
         fqsy = -y(PPICLF_JVY+j)/ppiclf_rprop(PPICLF_R_JTAUP,i)

         ! Gravity
         fbx  = 0.0
         fby  = -9.8

         ! set ydot for all PPICLF_LRS number of equations
         ydot(PPICLF_JX +j) = y(PPICLF_JVX +j)
         ydot(PPICLF_JY +j) = y(PPICLF_JVY +j)
         ydot(PPICLF_JVX+j) = fqsx+fbx
         ydot(PPICLF_JVY+j) = fqsy+fby
      enddo 
! evaluate ydot

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

      return
      end
!-----------------------------------------------------------------------
