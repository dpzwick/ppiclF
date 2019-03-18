!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(time_,y,ydot)
c
c     Set ydot(*) for system of d/dt y(*) = ydot(*)
c
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
c
c     Map Lagrangian quantity to Eulerian field:
c
c          map(*) = function of y(*), ydot(*), ydotc(*), rprop(*)
c
      real map(*)
      real y(*)
      real ydot(*)
      real ydotc(*)
      real rprop(*)

      return
      end
!-----------------------------------------------------------------------
