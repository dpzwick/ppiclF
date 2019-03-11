!-----------------------------------------------------------------------
      subroutine lpm_fun(time_,y,ydot)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      real time_
      real y(*)
      real ydot(*)

c setup interpolation
      call lpm_interpolate_setup
c setup interpolation

C interpolate fields
c     call lpm_interpolate_fld(LPM_R_JUX  , vx_e    )
c     call lpm_interpolate_fld(LPM_R_JUY  , vy_e    )
c     call lpm_interpolate_fld(LPM_R_JUZ  , vz_e    )
C interpolate fields

c evaluate ydot
      do i=1,lpm_npart
         ! striding solution y vector
         j = LPM_LRS*(i-1)

         ! fluid viscosity
         rmu   = 1.8E-5

         ! particle mass
         rmass = lpm_rprop(LPM_R_JVOLP,i)*lpm_rprop(LPM_R_JRHOP,i)

         ! Stokes drag force
         rdum  = 18.0*rmu/lpm_rprop(LPM_R_JDP,i)**2
         rdum  = rdum*lpm_rprop(LPM_R_JVOLP,i)
         fqsx  = rdum*(lpm_rprop(LPM_R_JUX,i) - y(LPM_JVX+j))
         fqsy  = rdum*(lpm_rprop(LPM_R_JUY,i) - y(LPM_JVY+j))
         fqsz  = rdum*(lpm_rprop(LPM_R_JUZ,i) - y(LPM_JVZ+j))

         ! Gravity
         fbx  = 0.0
         fby  = -9.8*rmass
         fbz  = 0.0

         ! set ydot for all LPM_SLN number of equations
         ydot(LPM_JX +j) = y(LPM_JVX +j)
         ydot(LPM_JY +j) = y(LPM_JVY +j)
         ydot(LPM_JZ +j) = y(LPM_JVZ +j)
         ydot(LPM_JVX+j) = (fqsx+fbx)/rmass
         ydot(LPM_JVY+j) = (fqsy+fby)/rmass
         ydot(LPM_JVZ+j) = (fqsz+fbz)/rmass
      enddo 
c evaluate ydot

c project fields
      call lpm_project
c project fields

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_project_map(map,y,ydot,ydotc,rprop)
c
c     map Lagrangian quantity to Eulerian field
c
      real map(*)
      real y(*)
      real ydot(*)
      real ydotc(*)
      real rprop(*)

      map(LPM_P_JPHIP) = rprop(LPM_R_JVOLP)   ! particle volume

      return
      end
!-----------------------------------------------------------------------
