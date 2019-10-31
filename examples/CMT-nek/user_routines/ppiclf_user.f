!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi, rmu, rhof, rmass, vmag, dp, rep, rphip, rphif, asndf,
     >       rmacr, rmachp, dvxdt, dvydt, rpr, rkappa, rcp_fluid,
     >       rcp_part, rmass_therm, rmass_add, rmach_rat
      real*8 fqsx, fqsy, fbx, fby, fcx, fcy, famx, famy, fdpdx, fdpdy,
     >       qq
      real*8 rcd1, rcd_mcr, rcd_std, rcd_M1, rcd_M2, rcd_am,
     >       beta, factor, f1M, f2M, f3M, C1, C2, C3, lrep

      integer*4 i
!
      rpi        = 4.0*atan(1.0)
      rmu        = 1.8E-5
      rcp_part   = 840.0
      rpr        = 0.72
      rcp_fluid  = 1006.0
      rkappa     = rcp_fluid*rmu/rpr

! evaluate ydot
      do i=1,ppiclf_npart
         ! Useful values
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)
         vmag  = (ppiclf_rprop(PPICLF_R_JUX,i)
     >           -ppiclf_y(PPICLF_JVX,i))**2
     >          +(ppiclf_rprop(PPICLF_R_JUY,i)
     >           -ppiclf_y(PPICLF_JVY,i))**2
         vmag  = sqrt(vmag)
         rhof  = ppiclf_rprop(PPICLF_R_JRHOF,i)
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         rep   = vmag*dp*rhof/rmu
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphip = min(0.3,rphip)
         rphif = 1.0-ppiclf_rprop(PPICLF_R_JPHIP,i)
         asndf = ppiclf_rprop(PPICLF_R_JCS,i)
         rmachp= vmag/asndf

         ! Quasi-steady force (Re_p and Ma_p corrections):
         !   Improved Drag Correlation for Spheres and Application 
         !   to Shock-Tube Experiments 
         !   - Parmar et al. (2010)
         !   - AIAA Journal
         if(rep .lt. 1E-12) then
            rcd1 = 1.0
         else 
            rmacr= 0.6 ! Critical rmachp no.
            rcd_mcr = (1.+0.15*rep**(0.684)) + 
     >                (rep/24.0)*(0.513/(1.+483./rep**(0.669)))
          if (rmachp .le. rmacr) then
             rcd_std = (1.+0.15*rep**(0.687)) + 
     >                (rep/24.0)*(0.42/(1.+42500./rep**(1.16)))
             rmach_rat = rmachp/rmacr
             rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
          else if (rmachp .le. 1.0) then
            rcd_M1 = (1.0+0.118*rep**0.813) +
     >                (rep/24.0)*0.69/(1.0+3550.0/rep**.793)
            C1 =  6.48
            C2 =  9.28
            C3 = 12.21
            f1M = -1.884 +8.422*rmachp -13.70*rmachp**2 +8.162*rmachp**3
            f2M = -2.228 +10.35*rmachp -16.96*rmachp**2 +9.840*rmachp**3
            f3M =  4.362 -16.91*rmachp +19.84*rmachp**2 -6.296*rmachp**3
            lrep = log(rep)
            factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
            rcd1 = rcd_mcr + (rcd_M1-rcd_mcr)*factor
          else if (rmachp .lt. 1.75) then
            rcd_M1 = (1.0+0.118*rep**0.813) +
     >              (rep/24.0)*0.69/(1.0+3550.0/rep**.793)
            rcd_M2 = (1.0+0.107*rep**0.867) +
     >              (rep/24.0)*0.646/(1.0+861.0/rep**.634)
            C1 =  6.48
            C2 =  8.93
            C3 = 12.21
            f1M = -2.963 +4.392*rmachp -1.169*rmachp**2 -0.027*rmachp**3
     >             -0.233*exp((1.0-rmachp)/0.011)
            f2M = -6.617 +12.11*rmachp -6.501*rmachp**2 +1.182*rmachp**3
     >             -0.174*exp((1.0-rmachp)/0.010)
            f3M = -5.866 +11.57*rmachp -6.665*rmachp**2 +1.312*rmachp**3
     >             -0.350*exp((1.0-rmachp)/0.012)
            lrep = log(rep)
            factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
            rcd1 = rcd_M1 + (rcd_M2-rcd_M1)*factor
          else
            rcd1 = (1.0+0.107*rep**0.867) +
     >                (rep/24.0)*0.646/(1.0+861.0/rep**.634)
          end if ! rmachp
         endif    ! rep
         beta = rcd1*3.0*rpi*rmu*dp
         ! Quasi-steady force (phi corrections):
         !   The Added Mass, Basset, and Viscous Drag Coefficients 
         !   in Nondilute Bubbly Liquids Undergoing Small-Amplitude 
         !   Oscillatory Motion
         !   - Sangani et al. (1991)
         !   - Phys. Fluids A
         beta = beta*(1.+2.*rphip)/(1-rphip)**3

         fqsx = beta
     >          *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))
         fqsy = beta
     >          *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))

         ! Added mass force (Ma_p corrections):
         !   On the unsteady inviscid force on cylinders and spheres 
         !   in subcritical compressible flow
         !   - Parmar et al. (2008)
         !   - Philos. Trans. R. Soc. London, Ser. A
         rcd_am = 1.0 + 1.8*rmachp**2 + 7.6*rmachp**4
         ! Added mass force (phi corrections):
         !   On the dispersed two-phase flow in the laminar flow regime
         !   - Zuber (1964)
         !   - Chem. Engng Sci.
         rcd_am = rcd_am*(1.0+2.0*rphip)
         rmass_add = rhof*ppiclf_rprop(PPICLF_R_JVOLP,i)*rcd_am

         famx = -rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDX,i)
         famy = -rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDY,i)

         dvxdt = ppiclf_ydot(PPICLF_JVX,i)
         dvydt = ppiclf_ydot(PPICLF_JVY,i)

         ! Pressure gradient force
         fdpdx = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDX,i)
         fdpdy = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >            ppiclf_rprop(PPICLF_R_JDPDY,i)
         
         ! Collision force:
         !  A discrete numerical model for granular assemblies
         !  - Cundall and Strack (1979)
         !  - Geotechnique
         call ppiclf_solve_NearestNeighbor(i)

         fcx  = ppiclf_ydotc(PPICLF_JVX,i)
         fcy  = ppiclf_ydotc(PPICLF_JVY,i)

         ! Quasi-steady heat transfer:
         !   Interaction of a planar shock wave with a dense particle 
         !   curtain: Modeling and experiments
         !   - Ling et al. (2012)
         !   - Physics of Fluids
         rmass_therm = rmass*rcp_part
         qq = 2.0*rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                           ppiclf_y(PPICLF_JT,i) )
         qq = qq*(1.0+0.3*rep**(0.5)*rpr**(1.0/3.0))

         ! set ydot for all PPICLF_SLN number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+famx+fdpdx+fcx)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+famy+fdpdy+fcy)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JT,i)  = qq/rmass_therm

!        ! Account for dv/dt added mass force for projection
         famx = famx - rmass_add*dvxdt
         famy = famy - rmass_add*dvydt

         ! Project hydrodynamic forces 
         ppiclf_ydotc(PPICLF_JVX,i)     = -(fqsx+famx)
         ppiclf_ydotc(PPICLF_JVY,i)     = -(fqsy+famy)
         ! Project work done by hydrodynamic forces:
         !   Inter-phase heat transfer and energy coupling in turbulent 
         !   dispersed multiphase flows
         !   - Ling et al. (2016)
         !   - Physics of Fluids
         ppiclf_rprop(PPICLF_R_JE,i) =fqsx*ppiclf_y(PPICLF_JVX,i)+
     >                                fqsy*ppiclf_y(PPICLF_JVY,i)+
     >                                famx*ppiclf_rprop(PPICLF_R_JUX,i)+
     >                                famy*ppiclf_rprop(PPICLF_R_JUY,i)
         ! Project heat transfer
         ppiclf_rprop(PPICLF_R_JE,i) = qq  +ppiclf_rprop(PPICLF_R_JE,i)
         ppiclf_rprop(PPICLF_R_JE,i) = -1.0*ppiclf_rprop(PPICLF_R_JE,i)
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
      dp_norm = 1.0/rprop(PPICLF_R_JDP)
      map(PPICLF_P_JPHIP)  = dp_norm*rprop(PPICLF_R_JVOLP)
      map(PPICLF_P_JFX)    = dp_norm*ydotc(PPICLF_JVX)
      map(PPICLF_P_JFY)    = dp_norm*ydotc(PPICLF_JVY)
      map(PPICLF_P_JE)     = dp_norm*rprop(PPICLF_R_JE)
      map(PPICLF_P_JPHIPU) = dp_norm*rprop(PPICLF_R_JVOLP)
     >                              *y(PPICLF_JVX)
      map(PPICLF_P_JPHIPV) = dp_norm*rprop(PPICLF_R_JVOLP)
     >                              *y(PPICLF_JVY)
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
!        rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
!    >          +rzdiff**2
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
!        rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x
     >            + (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y
!    >            + (yj(PPICLF_JVZ)-yi(PPICLF_JVZ))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
!        ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
!    >                              + rnmag*rn_12z

      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp

         ! give a bit larger collision threshold for walls
         rextra   = 0.0d0
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
!        rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
!    >          +rzdiff**2
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
!    >              -yi(PPICLF_JVZ)*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
!        ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
!    >                              + rnmag*rn_12z

      endif

      return
      end
!-----------------------------------------------------------------------
