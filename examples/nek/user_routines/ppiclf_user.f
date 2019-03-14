!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(istep,dt,time,y,ydot)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer istep
      real    dt
      real    time
      real    y(*)
      real    ydot(*)

c setup interpolation
      call ppiclf_solve_SetupInterp(istep,dt,time)
c setup interpolation

C interpolate fields
c     call ppiclf_solve_InterpField(PPICLF_R_JUX  , vx_e    )
c     call ppiclf_solve_InterpField(PPICLF_R_JUY  , vy_e    )
c     call ppiclf_solve_InterpField(PPICLF_R_JUZ  , vz_e    )
C interpolate fields

c evaluate ydot
      do i=1,ppiclf_npart
         ! striding solution y vector
         j = PPICLF_LRS*(i-1)

         ! fluid viscosity
         rmu   = 1.8E-5

         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Stokes drag force
         rdum  = 18.0*rmu/ppiclf_rprop(PPICLF_R_JDP,i)**2
         rdum  = rdum*ppiclf_rprop(PPICLF_R_JVOLP,i)
         fqsx  = rdum*(ppiclf_rprop(PPICLF_R_JUX,i) - y(PPICLF_JVX+j))
         fqsy  = rdum*(ppiclf_rprop(PPICLF_R_JUY,i) - y(PPICLF_JVY+j))
         fqsz  = rdum*(ppiclf_rprop(PPICLF_R_JUZ,i) - y(PPICLF_JVZ+j))

         ! Gravity
         fbx  = 0.0
         fby  = -9.8*rmass
         fbz  = 0.0

         ! set ydot for all PPICLF_SLN number of equations
         ydot(PPICLF_JX +j) = y(PPICLF_JVX +j)
         ydot(PPICLF_JY +j) = y(PPICLF_JVY +j)
         ydot(PPICLF_JZ +j) = y(PPICLF_JVZ +j)
         ydot(PPICLF_JVX+j) = (fqsx+fbx)/rmass
         ydot(PPICLF_JVY+j) = (fqsy+fby)/rmass
         ydot(PPICLF_JVZ+j) = (fqsz+fbz)/rmass
      enddo 
c evaluate ydot

c project fields
      call ppiclf_solve_ProjectParticleGrid
c project fields

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

      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)   ! particle volume

      return
      end
!-----------------------------------------------------------------------
