!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
#else
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4  npart
      real*8     y(*)
      real*8     rprop(*)
!
! Internal:
!
      integer*4 ppiclf_iglsum,ntotal
      external ppiclf_iglsum
!

      call ppiclf_prints('   *Begin AddParticles$')

      if (ppiclf_npart+npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   call ppiclf_exittr('Invalid number of particles$',
     >                      0.0,ppiclf_npart+npart)

      call ppiclf_printsi('      -Begin copy particles$',npart)

      ! First, append arrays onto existing arrays
      call ppiclf_copy(ppiclf_y(1,ppiclf_npart+1),
     >                 y,
     >                 npart*PPICLF_LRS)
      call ppiclf_copy(ppiclf_rprop(1,ppiclf_npart+1),
     >                 rprop,
     >                 npart*PPICLF_LRP)
      ppiclf_npart = ppiclf_npart + npart

      call ppiclf_printsi('      -Begin copy particles$',ppiclf_npart)

      if (.not. PPICLF_RESTART) then
         call ppiclf_prints('      -Begin ParticleTag$')
            call ppiclf_solve_SetParticleTag(npart)
         call ppiclf_prints('       End ParticleTag$')
      endif

      if (ppiclf_iglsum(ppiclf_npart,1).gt.0) then
         call ppiclf_prints('      -Begin CreateBin$')
            call ppiclf_comm_CreateBin
         call ppiclf_prints('       End CreateBin$')

         call ppiclf_prints('      -Begin FindParticle$')
            call ppiclf_comm_FindParticle
         call ppiclf_prints('       End FindParticle$')

         call ppiclf_prints('      -Begin MoveParticle$')
            call ppiclf_comm_MoveParticle
         call ppiclf_prints('       End MoveParticle$')

      endif

      call ppiclf_prints('    End AddParticles$')

      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop)
     > bind(C, name="ppiclc_solve_InitParticle")
#else
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
      integer*4  npart
      real*8     y(*)
      real*8     rprop(*)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitParticle$',0.0d0
     >   ,ppiclf_nid)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter must be before InitParticle$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)

      call ppiclf_prints('*Begin InitParticle$')

         call ppiclf_prints('   *Begin InitParam$')
            call ppiclf_solve_InitParam(imethod,ndim,iendian)
         call ppiclf_prints('    End InitParam$')
         
         if (.not. PPICLF_RESTART) then
            call ppiclf_solve_InitZero
            call ppiclf_solve_AddParticles(npart,y,rprop)
         endif

      call ppiclf_prints('   *Begin WriteParticleVTU$')
         call ppiclf_io_WriteParticleVTU('')
      call ppiclf_prints('    End WriteParticleVTU$')

!     call ppiclf_prints('   *Begin WriteBinVTU$')
!        call ppiclf_io_WriteBinVTU('')
!     call ppiclf_prints('    End WriteBinVTU$')

      call ppiclf_prints(' End InitParticle$')

            call ppiclf_io_OutputDiagGen

      PPICLF_LINIT = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParam(imethod,ndim,iendian)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
!
      if (imethod .eq. 0 .or. abs(imethod) .ge. 2)
     >   call ppiclf_exittr('Invalid integration method$',0.0d0,imethod)
      if (ndim .le. 1 .or. ndim .ge. 4)
     >   call ppiclf_exittr('Invalid problem dimension$',0.0d0,ndim)
      if (iendian .lt. 0 .or. iendian .gt. 1)
     >   call ppiclf_exittr('Invalid Endian$',0.0d0,iendian)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian

      ppiclf_filter       = -1   ! filt default for no projection

      ppiclf_iperiodic(1) = 1    ! periodic in x (== 0) 
      ppiclf_iperiodic(2) = 1    ! periodic in y (==0)
      ppiclf_iperiodic(3) = 1    ! periodic in z (==0)

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0d0
      ppiclf_time   = 0.0d0

      ppiclf_overlap    = .false.
      ppiclf_linit      = .false.
      ppiclf_lfilt      = .false.
      ppiclf_lfiltgauss = .false.
      ppiclf_lfiltbox   = .false.
      ppiclf_lintp      = .false.
      ppiclf_lproj      = .false.
      ppiclf_lsubbin    = .false.
      ppiclf_lsubsubbin = .false.
      if (PPICLF_INTERP .eq. 1)  ppiclf_lintp = .true.
      if (PPICLF_PROJECT .eq. 1) ppiclf_lproj = .true.

      ppiclf_xdrange(1,1) = -1E20
      ppiclf_xdrange(2,1) =  1E20
      ppiclf_xdrange(1,2) = -1E20
      ppiclf_xdrange(2,2) =  1E20
      ppiclf_xdrange(1,3) = -1E20
      ppiclf_xdrange(2,3) =  1E20

      ppiclf_d2chk(1) = 0.0d0
      ppiclf_d2chk(2) = 0.0d0
      ppiclf_d2chk(3) = 0.0d0

      ppiclf_n_bins(1) = 1
      ppiclf_n_bins(2) = 1
      ppiclf_n_bins(3) = 1

      ppiclf_bins_set(1) = 0
      ppiclf_bins_set(2) = 0
      ppiclf_bins_set(3) = 0

      ppiclf_bins_balance(1) = 0
      ppiclf_bins_balance(2) = 0
      ppiclf_bins_balance(3) = 0

      ppiclf_nwall    = 0
      ppiclf_iwallm   = 0

      PPICLF_INT_ICNT = 0

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
     > bind(C, name="ppiclc_solve_InitNeighborBin")
#else
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rwidth
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitNeighborBin$',0.0d0
     >                  ,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitNeighborBin$'
     >                  ,0.0d0,0)

      ppiclf_lsubsubbin = .true.

      ppiclf_d2chk(3) = rwidth

!BD:Testing passsing to picle
!      write(*,*) "Neigbor",rwidth     

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitTargetBins(str,n,balance)
     > bind(C, name="ppiclc_solve_InitTargetBins")
#else
      subroutine ppiclf_solve_InitTargetBins(str,n,balance)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 str
      integer*4 n
      integer*4 balance
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitTargetBins$'
     >                   ,0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitTargetBins$'
     >                  ,0.0d0,0)

      if (str == 'x' .or. str == 'X') then 
         ppiclf_n_bins(1) = n
         if (n .gt. 1) ppiclf_bins_set(1) = 1
         ppiclf_bins_balance(1) = balance
      elseif (str == 'y' .or. str == 'Y') then 
         ppiclf_n_bins(2) = n
         if (n .gt. 1) ppiclf_bins_set(2) = 1
         ppiclf_bins_balance(2) = balance
      elseif (str == 'z' .or. str == 'Z') then 
        if (ppiclf_ndim .lt. 3)
     >   call ppiclf_exittr('Dim must be 3 to use InitTargetBins on z$'
     >                   ,0.,ppiclf_ndim)
         ppiclf_n_bins(3) = n
         if (n .gt. 1) ppiclf_bins_set(3) = 1
         ppiclf_bins_balance(3) = balance
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetNeighborBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
!
      do i=1,ppiclf_npart
         ppiclf_nb_r(1,i) = floor((ppiclf_cp_map(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(2,i) = floor((ppiclf_cp_map(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_r(3,i) = floor((ppiclf_cp_map(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      enddo

      do i=1,ppiclf_npart_gp
         ppiclf_nb_g(1,i) = floor((ppiclf_rprop_gp(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(2,i) = floor((ppiclf_rprop_gp(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_g(3,i) = floor((ppiclf_rprop_gp(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitZero
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,j,ic,k,ie
!
      ic = 0
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ic = ic + 1
         ppiclf_y    (j,i) = 0.0d0
         ppiclf_ydot (j,i) = 0.0d0
         ppiclf_ydotc(j,i) = 0.0d0
         ppiclf_y1   (ic ) = 0.0d0
      enddo
      do j=1,PPICLF_LRP
         ppiclf_rprop(j,i) = 0.0d0
      enddo
      do j=1,PPICLF_LRP2
         ppiclf_rprop2(j,i) = 0.0d0
      enddo
      do j=1,PPICLF_LIP
         ppiclf_iprop(j,i) = 0
      enddo
      enddo
      ppiclf_npart = 0

      do ie=1,PPICLF_LEE
      do ic=1,PPICLF_LRP_INT
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
        ppiclf_int_fld(i,j,k,ic,ie) = 0.0d0
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_NearestNeighbor(i)
!
      implicit none
!
      include "PPICLF"
! 
! Input:
! 
      integer*4 i
! 
! Internal: 
! 
      real*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP)
      real*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
     >       dist_total
      integer*4 i_iim, i_iip, i_jjm, i_jjp, i_kkm, i_kkp, j, j_ii, j_jj,
     >          j_kk, jp
      real*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
     >       ab_dot_ac, ab_mag, ac_mag, zdist2
      integer*4 istride, k, kmax, kp, kkp, kk
! 
      i_iim = ppiclf_nb_r(1,i) - 1
      i_iip = ppiclf_nb_r(1,i) + 1
      i_jjm = ppiclf_nb_r(2,i) - 1
      i_jjp = ppiclf_nb_r(2,i) + 1
      i_kkm = ppiclf_nb_r(3,i) - 1
      i_kkp = ppiclf_nb_r(3,i) + 1

      dist2 = ppiclf_d2chk(3)**2

      do j=1,ppiclf_npart
         if (j .eq. i) cycle

         j_ii = ppiclf_nb_r(1,j)
         j_jj = ppiclf_nb_r(2,j)
         j_kk = ppiclf_nb_r(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         endif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         endif
         if (dist_total .gt. dist2) cycle

         call ppiclf_user_EvalNearestNeighbor(i,j,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_cp_map(1,j)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,j))

      enddo

      do j=1,ppiclf_npart_gp
         j_ii = ppiclf_nb_g(1,j)
         j_jj = ppiclf_nb_g(2,j)
         j_kk = ppiclf_nb_g(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         endif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         endif
         if (dist_total .gt. dist2) cycle

         jp = -1*j
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))

      enddo

      istride = ppiclf_ndim
      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            area = ppiclf_wall_n(4,j)
            rpz1 = ppiclf_cp_map(3,i)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         endif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         endif


         a_sum = 0.0d0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            endif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            if (rdist .gt. 2.0d0*ppiclf_d2chk(3)) goto 1511

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0d0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
               A(3) = ydum(3)
               B(3) = rpz1
               C(3) = rpz2
               AB(3) = B(3) - A(3)
               AC(3) = C(3) - A(3)

               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            endif
            a_sum = a_sum + tri_area
         enddo

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         if (a_sum .gt. rthresh*area) cycle

         jp = 0
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1511 continue
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitWall(xp1,xp2,xp3)
!
      implicit none
!
      include "PPICLF"
! 
! Input:
! 
      real*8 xp1(*)
      real*8 xp2(*)
      real*8 xp3(*)
!
! Internal:
!
      real*8 rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       a_sum, theta, tri_area, 
     >       ab_dot_ac, ab_mag, ac_mag, rise, run, rmag, 
     >       rpx3, rpy3, rpz3
      integer*4 istride, k, kmax, kp, kkp, kk
      real*8 A(3),B(3),C(3),AB(3),AC(3)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.d0,0)

      ppiclf_nwall = ppiclf_nwall + 1 

      if (ppiclf_nwall .gt. PPICLF_LWALL)
     >call ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.d0,ppiclf_nwall)

      istride = ppiclf_ndim
      a_sum = 0.0d0
      kmax = 2
      if (ppiclf_ndim .eq. 3) kmax = 3

      if (ppiclf_ndim .eq. 3) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp1(3)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(5,ppiclf_nwall) = xp2(2)
         ppiclf_wall_c(6,ppiclf_nwall) = xp2(3)
         ppiclf_wall_c(7,ppiclf_nwall) = xp3(1)
         ppiclf_wall_c(8,ppiclf_nwall) = xp3(2)
         ppiclf_wall_c(9,ppiclf_nwall) = xp3(3)

!TEMP Dump statments BAD
!        write(*,*) "PPWALL ", xp1(1),xp1(2),xp1(3)
!        write(*,*) "PPWALL ", xp2(1),xp2(2),xp2(3)
!        write(*,*) "PPWALL ", xp3(1),xp3(2),xp3(3)

         A(1) = (xp1(1) + xp2(1) + xp3(1))/3.0d0
         A(2) = (xp1(2) + xp2(2) + xp3(2))/3.0d0
         A(3) = (xp1(3) + xp2(3) + xp3(3))/3.0d0
      elseif (ppiclf_ndim .eq. 2) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(2)

         A(1) = (xp1(1) + xp2(1))/2.0d0
         A(2) = (xp1(2) + xp2(2))/2.0d0
         A(3) = 0.0d0
      endif

      ! compoute area:
      do k=1,kmax 
         kp = k+1
         if (kp .gt. kmax) kp = kp-kmax ! cycle
         
         kk   = istride*(k-1)
         kkp  = istride*(kp-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(kkp+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kkp+2,ppiclf_nwall)
         rpz2 = 0.0d0

         B(1) = rpx1
         B(2) = rpy1
         B(3) = 0.0d0
        
         C(1) = rpx2
         C(2) = rpy2
         C(3) = 0.0d0
        
         AB(1) = B(1) - A(1)
         AB(2) = B(2) - A(2)
         AB(3) = 0.0d0
        
         AC(1) = C(1) - A(1)
         AC(2) = C(2) - A(2)
         AC(3) = 0.0d0

         if (ppiclf_ndim .eq. 3) then
             rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
             rpz2 = ppiclf_wall_c(kkp+3,ppiclf_nwall)
             B(3) = rpz1
             C(3) = rpz2
             AB(3) = B(3) - A(3)
             AC(3) = C(3) - A(3)
        
             AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
             AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
             theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
             tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
         elseif (ppiclf_ndim .eq. 2) then
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
             tri_area = AB_MAG
         endif
         a_sum = a_sum + tri_area
      enddo
      
      ppiclf_wall_n(ppiclf_ndim+1,ppiclf_nwall) = a_sum

      ! wall normal:
      if (ppiclf_ndim .eq. 2) then

         rise = xp2(2) - xp1(2)
         run  = xp2(1) - xp1(1)

         rmag = sqrt(rise**2 + run**2)
         rise = rise/rmag
         run  = run/rmag
         
         ppiclf_wall_n(1,ppiclf_nwall) = rise
         ppiclf_wall_n(2,ppiclf_nwall) = -run

      elseif (ppiclf_ndim .eq. 3) then

         k  = 1
         kk = istride*(k-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 2
         kk = istride*(k-1)
         rpx2 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz2 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 3
         kk = istride*(k-1)
         rpx3 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy3 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz3 = ppiclf_wall_c(kk+3,ppiclf_nwall)
    
         A(1) = rpx2 - rpx1
         A(2) = rpy2 - rpy1
         A(3) = rpz2 - rpz1

         B(1) = rpx3 - rpx2
         B(2) = rpy3 - rpy2
         B(3) = rpz3 - rpz2

         ppiclf_wall_n(1,ppiclf_nwall) = A(2)*B(3) - A(3)*B(2)
         ppiclf_wall_n(2,ppiclf_nwall) = A(3)*B(1) - A(1)*B(3)
         ppiclf_wall_n(3,ppiclf_nwall) = A(1)*B(2) - A(2)*B(1)

         rmag = sqrt(ppiclf_wall_n(1,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(2,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(3,ppiclf_nwall)**2)

         ppiclf_wall_n(1,ppiclf_nwall) = ppiclf_wall_n(1,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(2,ppiclf_nwall) = ppiclf_wall_n(2,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(3,ppiclf_nwall) = ppiclf_wall_n(3,ppiclf_nwall)
     >                                  /rmag

      endif

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
     > bind(C, name="ppiclc_solve_InitPeriodicX")
#else
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 xl
      real*8 xr
! 
      if (xl .ge. xr)
     >call ppiclf_exittr('PeriodicX must have xl < xr$',xl,0)

      ppiclf_iperiodic(1) = 0

      ppiclf_xdrange(1,1) = xl
      ppiclf_xdrange(2,1) = xr

      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
     > bind(C, name="ppiclc_solve_InitPeriodicY")
#else
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 yl
      real*8 yr
! 
      if (yl .ge. yr)
     >call ppiclf_exittr('PeriodicY must have yl < yr$',yl,0)

      ppiclf_iperiodic(2) = 0

      ppiclf_xdrange(1,2) = yl
      ppiclf_xdrange(2,2) = yr

      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
     > bind(C, name="ppiclc_solve_InitPeriodicZ")
#else
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 zl
      real*8 zr
      real*8 zzr
      real*8 zzl    
! 
       !write(*,*) "BIG Z 1" 
      !zl = 0.0d0
      zzl = 0.0  
       !write(*,*) "BIG Z L" 
      !zr = 0.0001d0   
       zzr = 0.000125!0.0002 !0.0006  
        !write(*,*) "BIG Z R"
      if (1.eq.2) then ! (zl .ge. zr) then
      !write(*,*) zl,zr
       !write(*,*) "BIG Z BAD IF"     
      call ppiclf_exittr('PeriodicZ must have zl < zr$',zl,0)
      endif  
      if (ppiclf_ndim .lt. 3)
     > call ppiclf_exittr('Cannot do PeriodicZ if not 3D$',zl,0)

      !write(*,*) "BIG Z 2",zzl,zzr  
      ppiclf_iperiodic(3) = 0
      !write(*,*) "BIG Z 3"  
      ppiclf_xdrange(1,3) = zl
      !write(*,*) "BIG Z 4"  
      ppiclf_xdrange(2,3) = zr
        
      !write(*,*) "Per Z calling Init"  
      call ppiclf_solve_InitSolve

      !write(*,*) "End of PZ"  
      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
     > bind(C, name="ppiclc_solve_InitGaussianFilter")
#else
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      real*8    alpha
      integer*4 iwallm
! 
! Internal: 
! 
      real*8 rsig
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)
      if (iwallm .lt. 0 .or. iwallm .gt. 1)
     >call ppiclf_exittr('0 or 1 must be used to specify filter mirror$'
     >                  ,0.0d0,iwallm)

      ppiclf_filter = filt
      ppiclf_alpha  = alpha 
      ppiclf_iwallm = iwallm

      rsig             = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_alpha))

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT      = .true.
      PPICLF_LFILTGAUSS = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
     > bind(C, name="ppiclc_solve_InitBoxFilter")
#else
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      integer*4 iwallm
      integer*4 sngl_elem
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                   ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                   ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)

c     filt = sqrt(1.5d0*filt**2/log(2.0d0) + 1.0d0)

      ppiclf_filter = filt
      ppiclf_iwallm = iwallm

      ppiclf_d2chk(2)  = filt/2.0d0

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT    = .true.
      PPICLF_LFILTBOX = .true.

      ! option to only use the current element (filter width will be 
      ! ignored)
      ! Note that this assumes the element volume is that of
      ! a cuboid... will need to get a better way for general
      ! hexahedral element eventually
      if ( sngl_elem == 1 ) PPICLF_SNGL_ELEM = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetParticleTag(npart)
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 npart
! 
! Internal: 
! 
      integer*4 i
!
      do i=ppiclf_npart-npart+1,ppiclf_npart
         ppiclf_iprop(5,i) = ppiclf_nid 
         ppiclf_iprop(6,i) = ppiclf_cycle
         ppiclf_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
#else
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 istep
      integer*4 iostep
      real*8    dt
      real*8    time
! 
! Internal:
!
      logical iout
!
      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      if (ppiclf_imethod .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(iout)
      if (ppiclf_imethod .eq. -1) 
     >   call ppiclf_solve_IntegrateRK3s(iout)

      ! output files
      if (ppiclf_iostep .gt.0)then
      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then

         ! already wrote initial conditions
         if (ppiclf_cycle .ne. 0) then
            call ppiclf_io_WriteParticleVTU('')
            call ppiclf_io_WriteBinVTU('')
         endif

         if (ppiclf_lsubbin)
     >      call ppiclf_io_WriteSubBinVTU('')
      endif

      ! Output diagnostics
      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then
         call ppiclf_io_OutputDiagAll
      endif
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
!
! Output:
!
      logical iout
!
      ! save stage 1 solution
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = ppiclf_y(i,1)
      enddo

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_solve_SetYdot

         ! rk3 integrate
         do i=1,ndum
            ndum = PPICLF_NPART*PPICLF_LRS
            ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
         enddo
      enddo

      iout = .true.

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3s(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
      integer*4 icalld
      save      icalld
      data      icalld /0/
!
! Output:
!
      logical iout
!
      icalld = icalld + 1


      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = mod(icalld,nstage)
      if (istage .eq. 0) istage = 3
      iout = .false.
      if (istage .eq. nstage) iout = .true.

!BD: CHanged RK3 to match Rocflu's
!Check negative on third component

!Zero out for first stage
      if (istage .eq. 1) then
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = 0.0d0 
      enddo
      endif

      ! evaluate ydot
      call ppiclf_solve_SetYdot

      ! rk3 integrate CHECK SIGNS
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y(i,1) =  -ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                  + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                  + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
      enddo
!Store Current stage RHS for next stage's use
        do i=1,ndum
         ppiclf_y1(i) = ppiclf_ydot(i,1)
      enddo

!WAARNING: Experimental fix to keep particles unsure where to place this
!          command. Either before or after the storing of the current 
!          storage
        call ppiclf_solve_RemoveParticle      
!End Experimental fix

!BD: Original CODE Follows CMT-nek rk3 setup
      ! save stage 1 solution
!      if (istage .eq. 1) then
!      ndum = PPICLF_NPART*PPICLF_LRS
!      do i=1,ndum
!         ppiclf_y1(i) = ppiclf_y(i,1)
!      enddo
!      endif

      ! evaluate ydot
!      call ppiclf_solve_SetYdot

      ! rk3 integrate
!      ndum = PPICLF_NPART*PPICLF_LRS
!      do i=1,ndum
!         ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
!     >                  + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
!     >                  + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
!      enddo
!BD: original code END
      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetYdot
!
      implicit none
!
      include "PPICLF"
! 
      call ppiclf_solve_InitSolve
      call ppiclf_user_SetYdot
      call ppiclf_solve_RemoveParticle

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSolve
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, j
!
!      write(*,*) "1 Initsolve Start"  
      call ppiclf_comm_CreateBin
!        write(*,*) "2 Createbin" 
      call ppiclf_comm_FindParticle
!        write(*,*) "3 FineParticle" 
      call ppiclf_comm_MoveParticle
!        write(*,*) "4 MoveParticle" 
      if (ppiclf_overlap) 
     >   call ppiclf_comm_MapOverlapMesh
!        write(*,*) "5 MapOverlapMesh" 
!      if (ppiclf_lintp .and. ppiclf_int_icnt .ne. 0) 
      if ((ppiclf_lintp .and. ppiclf_int_icnt .ne. 0) .or.
     >    (ppiclf_lproj .and. ppiclf_sngl_elem))     
     >   call ppiclf_solve_InterpParticleGrid
!        write(*,*) "6 InterpPartGrid" 
      call ppiclf_solve_RemoveParticle
!        write(*,*) "7 RemovePart" 
      if (ppiclf_lsubsubbin .or. ppiclf_lproj) then
         call ppiclf_comm_CreateGhost
!        write(*,*) "8 CreateGhost" 
         call ppiclf_comm_MoveGhost
!        write(*,*) "9 MoveGhost" 
      endif

      if (ppiclf_lproj .and. ppiclf_overlap) 
     >   call ppiclf_solve_ProjectParticleGrid
!        write(*,*) "10 ProjectPart" 
      if (ppiclf_lsubsubbin) 
     >   call ppiclf_solve_SetNeighborBin
!        write(*,*) "11 SetNeighborBin" 
      ! Zero 
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ppiclf_ydotc(j,i) = 0.0d0
      enddo
      enddo
!        write(*,*) "12 InitSolve Done" 
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 j
!
      call ppiclf_solve_InitInterp
      do j=1,PPICLF_INT_ICNT
         call ppiclf_solve_InterpField(j)
      enddo
      call ppiclf_solve_FinalizeInterp

      call ppiclf_solve_PostInterp

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpFieldUser(jp,infld)
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 jp
      real*8 infld(*)
!
! Internal:
!
      integer*4 n
!
      if (PPICLF_INTERP .eq. 0)
     >call ppiclf_exittr(
     >     'No specified interpolated fields, set PPICLF_LRP_INT$',0.0d0
     >                   ,0)

      PPICLF_INT_ICNT = PPICLF_INT_ICNT + 1

      if (PPICLF_INT_ICNT .gt. PPICLF_LRP_INT)
     >   call ppiclf_exittr('Interpolating too many fields$'
     >                     ,0.0d0,PPICLF_INT_ICNT)
      if (jp .le. 0 .or. jp .gt. PPICLF_LRP)
     >   call ppiclf_exittr('Invalid particle array interp. location$'
     >                     ,0.0d0,jp)

      ! set up interpolation map
      PPICLF_INT_MAP(PPICLF_INT_ICNT) = jp

      ! copy to infld internal storage
      n = PPICLF_NEE*PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      call ppiclf_copy(ppiclf_int_fldu(1,1,1,1,PPICLF_INT_ICNT)
     >                ,infld(1),n)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitInterp
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 n, ie, npt_max, np, ndum
      real*8 tol, bb_t
!
      if (.not.ppiclf_overlap)
     >call ppiclf_exittr('Cannot interpolate unless overlap grid$',0.0d0
     >                   ,0)
      if (.not.ppiclf_lintp) 
     >call ppiclf_exittr('To interpolate, set PPICLF_LRP_PRO to ~= 0$'
     >                   ,0.0d0,0)

      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,1)
     >                   ,ppiclf_xm1b(1,1,1,1,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,2)
     >                   ,ppiclf_xm1b(1,1,1,2,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,3)
     >                   ,ppiclf_xm1b(1,1,1,3,ie),n)
      enddo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = 1
c     ndum    = ppiclf_neltb*n
      ndum    = ppiclf_neltb+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(ppiclf_fp_hndl
     >                         ,ppiclf_comm_nid
     >                         ,np ! only 1 rank on this comm
     >                         ,ppiclf_ndim
     >                         ,ppiclf_xm1bi(1,1,1,1,1)
     >                         ,ppiclf_xm1bi(1,1,1,1,2)
     >                         ,ppiclf_xm1bi(1,1,1,1,3)
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_neltb
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)


      ! copy MapOverlapMesh mapping from prior to communicating map
      ppiclf_neltbbb = ppiclf_neltbb
      do ie=1,ppiclf_neltbbb
         call ppiclf_icopy(ppiclf_er_mapc(1,ie),ppiclf_er_maps(1,ie)
     >             ,PPICLF_LRMAX)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(j)
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 jp
!
! Internal:
!
      integer*4 n, ie, iee, j
!
      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltbbb
         iee = ppiclf_er_mapc(1,ie)
         call ppiclf_copy(ppiclf_int_fld (1,1,1,j  ,ie)
     >                   ,ppiclf_int_fldu(1,1,1,iee,j ),n)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_FinalizeInterp
!
      implicit none
!
      include "PPICLF"
!
! Internal: 
!
      real*8 FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
      integer*4 nkey(2), nl, nii, njj, nxyz, nrr, ix, iy, iz, i, jp, ie
      logical partl
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltbbb
     >      ,PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld
     >      ,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltbbb
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld,nrr,nkey,2)

      ! find which cell particle is in locally
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim .eq. 3)
     >iz = 3

      call pfgslib_findpts(PPICLF_FP_HNDL           !   call pfgslib_findpts( ihndl,
     >        , ppiclf_iprop (1 ,1),PPICLF_LIP        !   $             rcode,1,
     >        , ppiclf_iprop (3 ,1),PPICLF_LIP        !   &             proc,1,
     >        , ppiclf_iprop (2 ,1),PPICLF_LIP        !   &             elid,1,
     >        , ppiclf_rprop2(1 ,1),PPICLF_LRP2       !   &             rst,ndim,
     >        , ppiclf_rprop2(4 ,1),PPICLF_LRP2       !   &             dist,1,
     >        , ppiclf_y     (ix,1),PPICLF_LRS        !   &             pts(    1),1,
     >        , ppiclf_y     (iy,1),PPICLF_LRS        !   &             pts(  n+1),1,
     >        , ppiclf_y     (iz,1),PPICLF_LRS ,PPICLF_NPART) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_INT_ICNT
         jp = PPICLF_INT_MAP(i)

         do ie=1,ppiclf_neltbbb
            call ppiclf_copy(fld(1,1,1,ie)
     >                      ,ppiclf_int_fld(1,1,1,i,ie),nxyz)
         enddo

         ! interpolate field locally
         call pfgslib_findpts_eval_local( PPICLF_FP_HNDL
     >                                  ,ppiclf_rprop (jp,1)
     >                                  ,PPICLF_LRP
     >                                  ,ppiclf_iprop (2,1)
     >                                  ,PPICLF_LIP
     >                                  ,ppiclf_rprop2(1,1)
     >                                  ,PPICLF_LRP2
     >                                  ,PPICLF_NPART
     >                                  ,fld)

      enddo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(PPICLF_FP_HNDL)

      ! Set interpolated fields to zero again
      PPICLF_INT_ICNT = 0

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_PostInterp
!
      implicit none
!
      include "PPICLF"
!
! Internal: 
!
      integer*4 rstride, istride
      parameter(rstride = 7 + PPICLF_LRP_INT)
      parameter(istride = 3)
      real*8 coord(rstride,PPICLF_LPART)
      integer*4 flag(istride,PPICLF_LPART)
      integer*4 fp_handle, i, j, k, npart
      external ppiclf_iglsum
      integer*4 ppiclf_iglsum
      integer*4 npt_max, np, ndum
      real*8 tol, bb_t
      integer*4 copy_back, jp, nxyz, ie
      real*8 xgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,ygrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,zgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
!
      ! Copy not found particles
      npart = 0
      do i=1,ppiclf_npart
         if (ppiclf_iprop(1,i) .eq. 2) then
            npart = npart + 1
            do j=1,ppiclf_ndim
               coord(j,npart) = ppiclf_y(j,i)
            enddo
         endif
      enddo

      if (ppiclf_iglsum(npart,1) .eq. 0) then
         return
      endif

      ! Copy grid indexing 
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_nee
      do i=1,nxyz
         xgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,1,ie)
         ygrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,2,ie)
         zgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,3,ie)
      enddo
      enddo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = ppiclf_np
c     ndum    = ppiclf_nee*n
      ndum    = ppiclf_nee+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(fp_handle
     >                         ,ppiclf_comm
     >                         ,np 
     >                         ,ppiclf_ndim
     >                         ,xgrid
     >                         ,ygrid
     >                         ,zgrid
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_nee
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)

      call pfgslib_findpts(fp_handle           !   call pfgslib_findpts( ihndl,
     >        , flag (1 ,1),istride        !   $             rcode,1,
     >        , flag (3 ,1),istride        !   &             proc,1,
     >        , flag (2 ,1),istride        !   &             elid,1,
     >        , coord(4 ,1),rstride       !   &             rst,ndim,
     >        , coord(7 ,1),rstride       !   &             dist,1,
     >        , coord(1,1) ,rstride        !   &             pts(    1),1,
     >        , coord(2,1) ,rstride        !   &             pts(  n+1),1,
     >        , coord(3,1) ,rstride ,npart) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_LRP_INT

         ! interpolate field (non-local)
         call pfgslib_findpts_eval( fp_handle
     >                                  ,coord (7+i,1)
     >                                  ,rstride
     >                                  ,flag (1,1)
     >                                  ,istride
     >                                  ,flag (3,1)
     >                                  ,istride
     >                                  ,flag (2,1)
     >                                  ,istride
     >                                  ,coord(4,1)
     >                                  ,rstride
     >                                  ,npart
     >                                  ,ppiclf_int_fldu(1,1,1,1,i))

      enddo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(fp_handle)

      ! Now copy particles back (assumes same ordering)
      k = 0
      do i=1,ppiclf_npart
         copy_back = 0
         if (ppiclf_iprop(1,i) .eq. 2) then
            k = k + 1
            if (flag(1,k) .lt. 2) then
               copy_back = 1
            endif
         endif

         if (copy_back .eq. 1) then
            ppiclf_iprop(1,i) = flag(1,k)
            do j=1,PPICLF_LRP_INT
               jp = PPICLF_INT_MAP(j)
               ppiclf_rprop(jp,i) = coord(7+j,k)
            enddo
         endif
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff(dt)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 dt
!


!BD:Rocflu's rk3 scheme

!Folowing form
        !rk3(1,:) = Temp storage i.e. Previous Stage RHS
        !rk3(2,:) = Temp storage i.e. Current Stage iteration
        !rk3(3,:) = Temp storage i.e. Current Stage RHS
      ppiclf_rk3coef(1,1) = 0.d00
      ppiclf_rk3coef(2,1) = 1.0d0
      ppiclf_rk3coef(3,1) = dt*8.0d0/15.0d0
      ppiclf_rk3coef(1,2) = dt*17.0d0/60.0d0
      ppiclf_rk3coef(2,2) = 1.0d0
      ppiclf_rk3coef(3,2) = dt*5.0d0/12.0d0
      ppiclf_rk3coef(1,3) = dt*5.0d0/12.0d0
      ppiclf_rk3coef(2,3) = 1.0d0
      ppiclf_rk3coef(3,3) = dt*3.0d0/4.0d0

!BD:Original Code This follows CMT-nek's rk 3 scheme
!      ppiclf_rk3coef(1,1) = 0.d00
!      ppiclf_rk3coef(2,1) = 1.0d0 
!      ppiclf_rk3coef(3,1) = dt
!      ppiclf_rk3coef(1,2) = 3.0d0/4.0d0
!      ppiclf_rk3coef(2,2) = 1.0d0/4.0d0 
!      ppiclf_rk3coef(3,2) = dt/4.0d0
!      ppiclf_rk3coef(1,3) = 1.0d0/3.0d0
!      ppiclf_rk3coef(2,3) = 2.0d0/3.0d0
!      ppiclf_rk3coef(3,3) = dt*2.0d0/3.0d0
!BD: Original Code END
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_MarkForRemoval(i)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
!
      ppiclf_iprop(1,i) = 3

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_RemoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >          iperiodicz,ndim, i, isl, isr, j, jchk, ic
!

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
      ndim       = ppiclf_ndim

      jj(1) = 1
      jj(2) = 2
      jj(3) = 3

      do i=1,ppiclf_npart
         isl = (i -1) * PPICLF_LRS + 1
         in_part(i) = 0
         if (ppiclf_iprop(1,i) .eq. 3) then
             write(*,*) "USER DP, should not be happening"   
            in_part(i) = -1 ! User removed particle
            goto 1513
         endif
         do j=0,ndim-1
            jchk = jj(j+1)
            if (ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
!     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
                  goto 1512
                endif
            endif
            if (ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) +
     &                     abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
!     &                     abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (ppiclf_iprop(1,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get herea
!              PCF dim xd1,xd2,val         
               write(*,*) "PCF", j+1,ppiclf_xdrange(1,j+1),
     >          ppiclf_xdrange(2,j+1),ppiclf_y(jchk,i)  
               
                  
            endif
 1512 continue
         enddo
 1513 continue
      enddo

      ic = 0
      do i=1,ppiclf_npart
!BAD DEBUG LOG for particles getting deleted
        if (in_part(i) .ne. 0) then
        write(*,*) "DelP in x y z, u v w", in_part(i), ppiclf_y(1,i),
     >   ppiclf_y(2,i) ,ppiclf_y(3,i),ppiclf_y(4,i),ppiclf_y(5,i),
     >   ppiclf_y(6,i)
        end if
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               isl = (i -1) * PPICLF_LRS + 1
               isr = (ic-1) * PPICLF_LRS + 1
               call ppiclf_copy
     >              (ppiclf_y     (1,ic),ppiclf_y(1,i)     ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_y1    (isr) ,ppiclf_y1(isl)    ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydot  (1,ic),ppiclf_ydot(1,i)  ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydotc (1,ic),ppiclf_ydotc(1,i) ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_rprop (1,ic),ppiclf_rprop(1,i) ,PPICLF_LRP)
               call ppiclf_copy
     >              (ppiclf_rprop2(1,ic),ppiclf_rprop2(1,i),PPICLF_LRP2)
               call ppiclf_icopy
     >              (ppiclf_iprop(1,ic) ,ppiclf_iprop(1,i) ,PPICLF_LIP)
            endif
         endif
      enddo
      if (ppiclf_npart .ne. ic) 
     >   write(*,*) "Tot DelP old, new",ppiclf_npart, ic  
      ppiclf_npart = ic

      return
      end
!----------------------------------------------------------------------
      subroutine ppiclf_solve_FindWallProject(rx2)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
       real*8 rx2(3)
!
! Internal:
!
      real*8 rnx,rny,rnz,rpx1,rpy1,rpz1,rpx2,rpy2,rpz2,rflip,rd,rdist
      integer*4 j, istride
!
      istride = ppiclf_ndim
      ppiclf_nwall_m = 0
      do j = 1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         rpx1 = rx2(1)
         rpy1 = rx2(2)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            rpz1 = rx2(3)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         endif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         endif

         rpx1 = ppiclf_wall_c(1,j)
         rpy1 = ppiclf_wall_c(2,j)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(istride+1,j)
         rpy2 = ppiclf_wall_c(istride+2,j)
         rpz2 = 0.0d0

         if (ppiclf_ndim .eq. 3) then
            rpz1 = ppiclf_wall_c(3,j)
            rpz2 = ppiclf_wall_c(istride+3,j)
         endif

         rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

         rdist = abs(rnx*rx2(1)+rny*rx2(2)
     >              +rnz*rx2(3)+rd)
         rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)
         rdist = rdist*2.0d0 ! Mirror

         if (rdist .gt. ppiclf_d2chk(2)) goto 1511

         ppiclf_nwall_m = ppiclf_nwall_m + 1

         ppiclf_xyz_mirror(1,ppiclf_nwall_m) = rx2(1) - rdist*rnx
         ppiclf_xyz_mirror(2,ppiclf_nwall_m) = rx2(2) - rdist*rny
         ppiclf_xyz_mirror(3,ppiclf_nwall_m) = 0.0d0
         if (ppiclf_ndim .eq. 3) then
            ppiclf_xyz_mirror(3,ppiclf_nwall_m) = rx2(3) - rdist*rnz
         endif

 1511 continue
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp
      logical partl, if3d
      integer*4 nkey(2), nxyz, nxyzdum, i, j, k, idum, ic, ip, iip, jjp,
     >          kkp, ilow, ihigh, jlow, jhigh, klow, khigh, ie, jj, j1,
     >          neltbc, ndum, nl, nii, njj, nrr, nlxyzep, iee, ndumdum
      real*8 pi, d2chk2_sq, rdum, multfci, rsig, rdist2, rexp, rx2(3),
     >       rx22, ry22, rz22, rtmp2, evol
!
      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      nxyzdum = nxyz*PPICLF_LRP_PRO*PPICLF_LEE
      do i=1,nxyzdum
         ppiclf_pro_fldb(i,1,1,1,1) = 0.0d0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d) ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d0)
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         if (ppiclf_sngl_elem) then
           multfci = 1.0d0
           rdum = multfci
         else
           multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
           if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
           rdum = multfci
         endif
      endif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d) 
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         do j=idum+1,idum+PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         enddo
                    
         iproj(1,ip)  = ppiclf_iprop(8,ip)
         iproj(2,ip)  = ppiclf_iprop(9,ip)
         if (if3d)
     >   iproj(3,ip)  = ppiclf_iprop(10,ip)
         iproj(4,ip)  = ppiclf_iprop(11,ip)
      enddo

      if (.not. ppiclf_sngl_elem) then
  
        ! ghost particles
        do ip=1,ppiclf_npart_gp
  
           rproj(1 ,ip+ppiclf_npart) = rdum
           rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
           rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
           if (if3d) 
     >     rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)
  
           idum = PPICLF_LRS+PPICLF_LRP
           ic = 4
           do j=idum+1,idum+PPICLF_LRP_GP
              ic = ic + 1
              rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
           enddo
                      
           iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
           iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
           if (if3d)
     >     iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
           iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
        enddo
  
        ndum = ppiclf_npart+ppiclf_npart_gp
  
        do ip=1,ndum
           iip      = iproj(1,ip)
           jjp      = iproj(2,ip)
           if (if3d)
     >     kkp      = iproj(3,ip)
           ndumdum  = iproj(4,ip)
  
           ilow  = iip-1
           ihigh = iip+1
           jlow  = jjp-1
           jhigh = jjp+1
           if (if3d) then
              klow  = kkp-1
              khigh = kkp+1
           endif
  
           ! Find if particle near wall and should mirror itself
           if (ppiclf_iwallm .eq. 1) then
              rx2(1) = rproj(2,ip)
              rx2(2) = rproj(3,ip)
              rx2(3) = rproj(4,ip)
              call ppiclf_solve_FindWallProject(rx2)
           endif
  
           do ie=1,ppiclf_neltb
  
                 if (ppiclf_el_map(1,ie) .gt. ndumdum) exit
                 if (ppiclf_el_map(2,ie) .lt. ndumdum) cycle 
           
                 if (ppiclf_el_map(3,ie) .gt. ihigh) cycle
                 if (ppiclf_el_map(4,ie) .lt. ilow)  cycle
                 if (ppiclf_el_map(5,ie) .gt. jhigh) cycle
                 if (ppiclf_el_map(6,ie) .lt. jlow)  cycle
                 if (if3d) then
                 if (ppiclf_el_map(7,ie) .gt. khigh) cycle
                 if (ppiclf_el_map(8,ie) .lt. klow)  cycle
                 endif
  
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              if (ppiclf_modgp(i,j,k,ie,4).ne.ndumdum) cycle
  
              rdist2  = (ppiclf_xm1b(i,j,k,1,ie) - rproj(2,ip))**2 +
     >                  (ppiclf_xm1b(i,j,k,2,ie) - rproj(3,ip))**2
              if(if3d) rdist2 = rdist2 +
     >                  (ppiclf_xm1b(i,j,k,3,ie) - rproj(4,ip))**2
  
              if (rdist2 .gt. d2chk2_sq) cycle
  
              rexp = 1.0d0
              if (ppiclf_lfiltgauss)
     >           rexp = exp(rdist2*rproj(1,ip))
  
              ! add wall effects
              if (ppiclf_iwallm .eq. 1) then
                 do jj=1,ppiclf_nwall_m
                    rx22 = (ppiclf_xm1b(i,j,k,1,ie) 
     >                     -ppiclf_xyz_mirror(1,jj))**2
                    ry22 = (ppiclf_xm1b(i,j,k,2,ie)
     >                     -ppiclf_xyz_mirror(2,jj))**2
                    rtmp2 = rx22 + ry22
                    if (if3d) then
                       rz22 = (ppiclf_xm1b(i,j,k,3,ie)
     >                        -ppiclf_xyz_mirror(3,jj))**2
                       rtmp2 = rtmp2 + rz22
                    endif
                    if (ppiclf_lfiltgauss) then
                       rexp = rexp + exp(rtmp2*rproj(1,ip))
                    else
                       rexp = rexp + 1.0d0
                    endif
                 enddo
              endif
  
              
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              enddo
           enddo
           enddo
           enddo
           enddo
        enddo
      ! sngl elem
      else

!BD: This is where rproj gets stored in fld, this is one of the steps
!that should be tracked
 

        do ip=1,ppiclf_npart
!         write(*,*) "BD:Proj; ip,ipTot,elem+1,ElemTOT",
!     >    ip,ppiclf_npart,ppiclf_iprop(2,ip)+1,ppiclf_neltb
        
           do ie=1,ppiclf_neltb
  
             ! Only use the current (single) element from findpts
             ! Note that this assumes the element volume is that of
             ! a cuboid... will need to get a better way for general
             ! hexahedral element eventually
             if (ie .ne. ppiclf_iprop(2,ip)+1) cycle
!BD: Check if paricles are being found for projection!
!             write(*,*) "Part found in sngl_elm for proj",ip,ie   

             evol = (ppiclf_xm1b(PPICLF_LEX,1,1,1,ie) 
     >             - ppiclf_xm1b(1,1,1,1,ie))
             evol = evol
     >            * (ppiclf_xm1b(1,PPICLF_LEY,1,2,ie) 
     >             - ppiclf_xm1b(1,1,1,2,ie))
             if (if3d)
     >       evol = evol
     >            * (ppiclf_xm1b(1,1,PPICLF_LEZ,3,ie) 
     >             - ppiclf_xm1b(1,1,1,3,ie))
             rexp = 1.0 / evol
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              enddo
           enddo
           enddo
           enddo
           enddo
        enddo
      endif

      ! now send xm1b to the processors in nek that hold xm1

      neltbc = ppiclf_neltb
      ndum = PPICLF_LRMAX*neltbc
      call ppiclf_icopy(ppiclf_er_mapc,ppiclf_er_map,ndum)
      do ie=1,neltbc
         ppiclf_er_mapc(5,ie) = ppiclf_er_mapc(2,ie)
         ppiclf_er_mapc(6,ie) = ppiclf_er_mapc(2,ie)
      enddo
      nl = 0
      nii = PPICLF_LRMAX
      njj = 6
      nrr = nxyz*PPICLF_LRP_PRO
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc,
     >   PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,neltbc
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,nkey,2)

      ! add the fields from the bins to ptw array
      nlxyzep = nxyz*PPICLF_LEE*PPICLF_LRP_PRO
      do i=1,nlxyzep
         PPICLF_PRO_FLD(i,1,1,1,1) = 0.0d0
      enddo


      do ie=1,neltbc
         iee = ppiclf_er_mapc(1,ie)
         do ip=1,PPICLF_LRP_PRO
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
           PPICLF_PRO_FLD(i,j,k,iee,ip) = PPICLF_PRO_FLD(i,j,k,iee,ip) +
     >                                    PPICLF_PRO_FLDB(i,j,k,ip,ie)

         enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp, nxyz, nxyzdum,
     >          idum, jdum, kdum, ic, i, j, k, ip, ndum, il, ir, jl, jr,
     >          kl, kr, jj, j1, iip, jjp, kkp
      logical if3d
      real*8 pi, d2chk2_sq, rdum, rsig, multfci, rexp, rdist2
!

      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0d0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d)
     >ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d0)
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
         if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
      endif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         do j=idum+1,idum+PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         enddo

         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip) = 
     >       floor( (rproj(4,ip) - ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ! ghost particles
      do ip=1,ppiclf_npart_gp

         rproj(1 ,ip+ppiclf_npart) = rdum
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         do j=idum+1,idum+PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
         enddo
                    
         iproj(1,ip+ppiclf_npart) = 
     >     floor((rproj(2,ip+ppiclf_npart)-ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip+ppiclf_npart) = 
     >     floor((rproj(3,ip+ppiclf_npart)-ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip+ppiclf_npart) = 
     >     floor((rproj(4,ip+ppiclf_npart)-ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp

      if (ppiclf_lfiltgauss) then
         idum = floor(ppiclf_filter/2.0d0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         jdum = floor(ppiclf_filter/2.0d0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         kdum = 999999999
         if (if3d)
     >   kdum = floor(ppiclf_filter/2.0d0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
      endif

      if (ppiclf_lfiltbox) then
         idum = ppiclf_ngrids/2+1
         jdum = ppiclf_ngrids/2+1
         kdum = 999999999
         if (if3d)
     >   kdum = ppiclf_ngrids/2+1
      endif

      do ip=1,ndum
         iip = iproj(1,ip)
         jjp = iproj(2,ip)
         if (if3d)
     >   kkp = iproj(3,ip)

         il  = max(1     ,iip-idum)
         ir  = min(ppiclf_bx,iip+idum)
         jl  = max(1     ,jjp-jdum)
         jr  = min(ppiclf_by,jjp+jdum)
         kl  = 1
         kr  = 1
         if (if3d) then
         kl  = max(1     ,kkp-kdum)
         kr  = min(ppiclf_bz,kkp+kdum)
         endif

c        do k=kl,kr
c        do j=jl,jr
c        do i=il,ir
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
            if(if3d) rdist2 = rdist2 +
     >                (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = 1.0d0
            if (ppiclf_lfiltgauss)
     >         rexp = exp(rdist2*rproj(1,ip))

            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_grid_fld(i,j,k,jj) = 
     >                         ppiclf_grid_fld(i,j,k,jj) 
     >                       + sngl(rproj(j1,ip)*rexp)
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_GetProFldIJKEF(i,j,k,e,m,fld)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i,j,k,e,m
!
! Output:
!
      real*8 fld
!
      fld = ppiclf_pro_fld(i,j,k,e,m)

      return
      end
!-----------------------------------------------------------------------
