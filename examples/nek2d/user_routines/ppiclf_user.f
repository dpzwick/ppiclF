!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(time,y,ydot)
#include "PPICLF"

      real    time
      real    y(*)
      real    ydot(*)

C interpolate fields
      call ppiclf_solve_InitInterp
         call ppiclf_solve_InterpField(PPICLF_R_JPHIP,ppiclf_pro_fld)
      call ppiclf_solve_FinalizeInterp
C interpolate fields

c evaluate ydot
      do i=1,ppiclf_npart
         ! striding solution y vector
         j = PPICLF_LRS*(i-1)

         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Gravity
         fbx  = 0.0
         fby  = -9.8*rmass

         call ppiclf_solve_NearestNeighbor(i)

         fcx  = ppiclf_ydotc(PPICLF_JVX,i)
         fcy  = ppiclf_ydotc(PPICLF_JVY,i)

         ! set ydot for all PPICLF_SLN number of equations
         ydot(PPICLF_JX +j) = y(PPICLF_JVX +j)
         ydot(PPICLF_JY +j) = y(PPICLF_JVY +j)
         ydot(PPICLF_JVX+j) = (fbx+fcx)/rmass
         ydot(PPICLF_JVY+j) = (fby+fcy)/rmass
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

      ! particle volume divided by particle diameter for 2d
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)/rprop(PPICLF_R_JDP)

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

      rksp  = 10000.0
      erest = 0.3
      
      rpi2  =  9.869604401089358

      if (j .ne. 0) then
         rthresh  = 0.5*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
         rmult = 1./sqrt(1./rm1+1./rm2)
         eta   = 2.*sqrt(rksp)*log(erest)/sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x +
     >           (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y 

         rv12_mage = rv12_mag*eta
         rksp_max = rksp*rdelta12
         rnmag = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
      elseif (j .eq. 0) then

         rthresh  = 0.5*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.*sqrt(rksp)*log(erest)/sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -1.0*(yi(PPICLF_JVX)*rn_12x +
     >                    yi(PPICLF_JVY)*rn_12y)

         rv12_mage = rv12_mag*eta
         rksp_max = rksp*rdelta12
         rnmag = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y

      endif

      return
      end
!-----------------------------------------------------------------------
