!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi, rmu, rhof, rmass, vmag, dp, rep, rphip, rphif, 
     >       beta, reps, cds, fqsx, fqsy, fbx, fby, fcx, fcy,
     >       beta1, beta2, rs, rdum, fdpdx, fdpdy, fqsz, fdpdz, fcz

      integer*4 i
!
      rpi  = 4.0*atan(1.0)
      rmu  = 1.8E-5

! evaluate ydot
      do i=1,ppiclf_npart
         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Gidaspow drag
         vmag  = (ppiclf_rprop(PPICLF_R_JUX,i)
     >           -ppiclf_y(PPICLF_JVX,i))**2
     >          +(ppiclf_rprop(PPICLF_R_JUY,i)
     >           -ppiclf_y(PPICLF_JVY,i))**2
     >          +(ppiclf_rprop(PPICLF_R_JUZ,i)
     >           -ppiclf_y(PPICLF_JVZ,i))**2
         vmag  = sqrt(vmag)
         rhof  = ppiclf_rprop(PPICLF_R_JRHOF,i)
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         rep   = vmag*dp*rhof/rmu
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphif = 1.0-ppiclf_rprop(PPICLF_R_JPHIP,i)

         if (rphip .gt. 0.2) then
            beta = 150.*rphip*rmu/rphif/dp**2 + 1.75*rhof*vmag/dp
         else
            reps = rep*rphif
            reps = max(reps,1E-12)
            if (reps .le. 1E3) then
               cds = 24.0/reps*(1.0+0.15*reps**(0.687))
            else
               cds = 0.44
            endif
            beta = 0.75*cds*rphif*rhof*vmag*rphif**(-2.65)/dp
         endif

         fqsx = beta*ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))
         fqsy = beta*ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))
         fqsz = beta*ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JUZ,i)-ppiclf_y(PPICLF_JVZ,i))

         fdpdx = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDX,i)
         fdpdy = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDY,i)
         fdpdz = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDZ,i)
         
         ! Collision search
         call ppiclf_solve_NearestNeighbor(i)

         fcx  = ppiclf_ydotc(PPICLF_JVX,i)
         fcy  = ppiclf_ydotc(PPICLF_JVY,i)
         fcz  = ppiclf_ydotc(PPICLF_JVZ,i)

         ! set ydot for all PPICLF_SLN number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+fdpdx+fcx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+fdpdy+fcy)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+fdpdz+fcz)/rmass

         ppiclf_ydotc(PPICLF_JVX,i) = -fqsx
         ppiclf_ydotc(PPICLF_JVY,i) = -fqsy
         ppiclf_ydotc(PPICLF_JVZ,i) = -fqsz
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
! Internal:
!
      real*8 dp_norm
!
      ! particle volume divided by particle diameter for 2d
!     dp_norm = 1./rprop(PPICLF_R_JDP)
!     map(PPICLF_P_JPHIP) = dp_norm*rprop(PPICLF_R_JVOLP)
!     map(PPICLF_P_JFX)   = dp_norm*ydotc(PPICLF_JVX)
!     map(PPICLF_P_JFY)   = dp_norm*ydotc(PPICLF_JVY)
      ! not in 3D though
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ)

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
! Internal:
!
      real*8 ksp,erest
      common /ucollision/ ksp,erest

      real*8 rpi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra
!
      rpi2  =  9.869604401089358d0

      ! other particles
      if (j .ne. 0) then
         rthresh  = 0.5d0*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
     >          +rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
         rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+rpi2)
     >           *rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x
     >            + (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y
     >            + (yj(PPICLF_JVZ)-yi(PPICLF_JVZ))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z

      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp

         ! give a bit larger collision threshold for walls
         rextra   = 0.0d0
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
     >          +rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.0d0*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -yi(PPICLF_JVX)*rn_12x
     >              -yi(PPICLF_JVY)*rn_12y
     >              -yi(PPICLF_JVZ)*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z

      endif

      return
      end
!-----------------------------------------------------------------------
