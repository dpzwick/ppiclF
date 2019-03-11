!-----------------------------------------------------------------------
      subroutine lpm_fun(time_,y,ydot)
c
c     Set ydot(*) for system of d/dt y(*) = ydot(*)
c
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real time_
      real y(*)
      real ydot(*)

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_project_map(map,y,ydot,ydotc,rprop)
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
