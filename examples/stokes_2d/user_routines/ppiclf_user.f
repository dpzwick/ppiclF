!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(time,y,ydot)
#include "PPICLF"

      real    time
      real    y(*)
      real    ydot(*)

c evaluate ydot
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
c evaluate ydot

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
c
c     map Lagrangian quantity to Eulerian field
c
      real map(*)
      real y(*)
      real ydot(*)
      real ydotc(*)
      real rprop(*)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
#include "PPICLF"
c
c     called for every j nearest neighbor of particle i
c
      integer i
      integer j
      real yi(*)     ! PPICLF_LRS
      real rpropi(*) ! PPICLF_LRP
      real yj(*)     ! PPICLF_LRS
      real rpropj(*) ! PPICLF_LRP

      return
      end
!-----------------------------------------------------------------------
