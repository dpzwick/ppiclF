!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot(time,y,ydot)
#include "PPICLF"

      real    time
      real    y(*)
      real    ydot(*)

      ! External from Nek5000
      real fld_to_interp(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE,
     >                   PPICLF_LRP_INT)
      common /interp_fld_nek/ fld_to_interp
      ! External from Nek5000

C interpolate fields
      call ppiclf_solve_InitInterp
         call ppiclf_solve_InterpField
     >                 (PPICLF_R_JPHIP,fld_to_interp(1,1,1,1,1))
         call ppiclf_solve_InterpField
     >                 (PPICLF_R_JUX,fld_to_interp(1,1,1,1,2)  )
         call ppiclf_solve_InterpField
     >                 (PPICLF_R_JUY,fld_to_interp(1,1,1,1,3)  )
      call ppiclf_solve_FinalizeInterp
C interpolate fields

      rpi  = 4.0*atan(1.0)
      rmu  = 1.8E-5
      rhof = 1.205

c evaluate ydot
      do i=1,ppiclf_npart
         ! striding solution y vector
         j = PPICLF_LRS*(i-1)

         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Gidaspow drag below
         vmag  = sqrt((ppiclf_rprop(PPICLF_R_JUX,i)-y(PPICLF_JVX+j))**2
     >               +(ppiclf_rprop(PPICLF_R_JUY,i)-y(PPICLF_JVY+j))**2)
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         rep   = vmag*dp*rhof/rmu
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphif = 1.0-ppiclf_rprop(PPICLF_R_JPHIP,i)

         if (rphip .gt. 0.2) then
            beta = 150.*rphip*rmu/rphif/dp**2 + 1.75*rhof*vmag/dp
         else
            reps = rep*rphif
            reps = max(reps,1E-1)
            if (reps .le. 1E3) then
               cds = 24.0/reps*(1.0+0.15*reps**(0.687))
            else
               cds = 0.44
            endif
            beta = 0.75*cds*rphif*rhof*vmag*rphif**(-2.65)/dp
         endif


         fqsx = beta*ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JUX,i)-y(PPICLF_JVX+j))
         fqsy = beta*ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JUY,i)-y(PPICLF_JVY+j))

         ! Gravity
         fbx  = 0.0
         fby  = -9.8*rmass

         call ppiclf_solve_NearestNeighbor(i)

         fcx  = ppiclf_ydotc(PPICLF_JVX,i)
         fcy  = ppiclf_ydotc(PPICLF_JVY,i)

         ! set ydot for all PPICLF_SLN number of equations
         ydot(PPICLF_JX +j) = y(PPICLF_JVX +j)
         ydot(PPICLF_JY +j) = y(PPICLF_JVY +j)
         ydot(PPICLF_JVX+j) = (fbx+fqsx+fcx)/rmass
         ydot(PPICLF_JVY+j) = (fby+fqsy+fcy)/rmass

         ! Store drag force for coupling with fluid
         ppiclf_ydotc(PPICLF_JVX,i) = fqsx
         ppiclf_ydotc(PPICLF_JVY,i) = fqsy
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
      dp_norm = 1./rprop(PPICLF_R_JDP)
      map(PPICLF_P_JPHIP) = dp_norm*rprop(PPICLF_R_JVOLP)
      map(PPICLF_P_JFX)   = dp_norm*ydotc(PPICLF_JVX)
      map(PPICLF_P_JFY)   = dp_norm*ydotc(PPICLF_JVY)

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

      ! For user implemented collision model
      real ksp,erest
      common /external_user_collsion/ ksp,erest
      ! For user implemented collision model
      
      rpi2  =  9.869604401089358

      ! other particles
      if (j .ne. 0) then
         rthresh  = 0.5*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
         rmult = 1./sqrt(1./rm1+1./rm2)
         eta   = 2.*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x +
     >           (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y

      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp

         ! give a bit larger collision threshold for walls
         rextra   = 0.0
         rthresh  = (0.5+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1./rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -1.0*(yi(PPICLF_JVX)*rn_12x +
     >                    yi(PPICLF_JVY)*rn_12y )

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         PPICLF_YDOTC(PPICLF_JVX,i) = PPICLF_YDOTC(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         PPICLF_YDOTC(PPICLF_JVY,i) = PPICLF_YDOTC(PPICLF_JVY,i)
     >                              + rnmag*rn_12y

      endif

      return
      end
!-----------------------------------------------------------------------
