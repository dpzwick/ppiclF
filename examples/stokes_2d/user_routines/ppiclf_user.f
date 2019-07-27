!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 fqsx,fqsy,fbx,fby
      integer*4 i
!
! evaluate ydot
      do i=1,ppiclf_npart
         ! Stokes drag
         fqsx = -ppiclf_y(PPICLF_JVX,i)/ppiclf_rprop(PPICLF_R_JTAUP,i)
         fqsy = -ppiclf_y(PPICLF_JVY,i)/ppiclf_rprop(PPICLF_R_JTAUP,i)

         ! Gravity
         fbx  = 0.0d0
         fby  = -9.8d0

         ! set ydot for all PPICLF_LRS number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JVX,i) = fqsx+fbx
         ppiclf_ydot(PPICLF_JVY,i) = fqsy+fby
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
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)    
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)    
      real*8 rpropj(PPICLF_LRP)
!

      return
      end
!-----------------------------------------------------------------------
