!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y)
#include "PPICLF"

      integer  imethod
      integer  ndim
      integer  iendian
      integer  npart
      real     y(*)

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitParticle$',0.
     >   ,ppiclf_nid)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter must be before InitParticle$',0.,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.,0)

      call ppiclf_prints('*Begin InitParticle$')

         call ppiclf_prints('   *Begin InitParam$')
            call ppiclf_solve_InitParam(imethod,ndim,iendian,npart)
         call ppiclf_prints('    End InitParam$')
         
            call ppiclf_copy(ppiclf_y,y,ppiclf_npart)
         
         call ppiclf_prints('   *Begin ParticleTag$')
            call ppiclf_solve_InitParticleTag
            call ppiclf_solve_SetParticleTag
         call ppiclf_prints('    End ParticleTag$')

         call ppiclf_prints('   *Begin CreateBin$')
            call ppiclf_comm_CreateBin
         call ppiclf_prints('    End CreateBin$')

         call ppiclf_prints('   *Begin FindParticle$')
            call ppiclf_comm_FindParticle
         call ppiclf_prints('    End FindParticle$')

         call ppiclf_prints('   *Begin MoveParticle$')
            call ppiclf_comm_MoveParticle
         call ppiclf_prints('    End MoveParticle$')

         call ppiclf_prints('   *Begin WriteParticleVTU$')
            call ppiclf_io_WriteParticleVTU('')
         call ppiclf_prints('    End WriteParticleVTU$')

         call ppiclf_prints('   *Begin WriteBinVTU$')
            call ppiclf_io_WriteBinVTU('')
         call ppiclf_prints('    End WriteBinVTU$')
         
      call ppiclf_prints(' End InitParticle$')

            call ppiclf_io_OutputDiagGen

      PPICLF_LINIT = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParam(imethod,ndim,iendian,npart)
#include "PPICLF"
      integer  imethod
      integer  ndim
      integer  iendian
      integer  npart

      if (imethod .le. 0 .or. imethod .ge. 3)
     >   call ppiclf_exittr('Invalid integration method$',0.0,imethod)
      if (ndim .le. 1 .or. ndim .ge. 4)
     >   call ppiclf_exittr('Invalid problem dimension$',0.0,ndim)
      if (iendian .lt. 0 .or. iendian .gt. 1)
     >   call ppiclf_exittr('Invalid Endian$',0.0,iendian)
      if (npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   call ppiclf_exittr('Invalid number of particles$',0.0,npart)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian
      ppiclf_npart        = npart

      ppiclf_filter       = -1   ! filt default for no projection

      ppiclf_iperiodic(1) = 1    ! periodic in x (== 0) 
      ppiclf_iperiodic(2) = 1    ! periodic in y (==0)
      ppiclf_iperiodic(3) = 1    ! periodic in z (==0)

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0
      ppiclf_time   = 0.0

      ppiclf_restart    = .false.
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

      ppiclf_xdrange(1,1) = -1E10
      ppiclf_xdrange(2,1) =  1E10
      ppiclf_xdrange(1,2) = -1E10
      ppiclf_xdrange(2,2) =  1E10
      ppiclf_xdrange(1,3) = -1E10
      ppiclf_xdrange(2,3) =  1E10

      ppiclf_d2chk(1) = 0
      ppiclf_d2chk(2) = 0
      ppiclf_d2chk(3) = 0

      ppiclf_nbin_dir(1) = 0
      ppiclf_nbin_dir(2) = 0
      ppiclf_nbin_dir(3) = 0

      ppiclf_nwall = 0

      PPICLF_INT_ICNT = -1

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
#include "PPICLF"

      real rwidth

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitNeighborBin$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitNeighborBin$'
     >                  ,0.,0)

      ppiclf_lsubsubbin = .true.

      ppiclf_d2chk(3) = rwidth

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSuggestedDir(str)
#include "PPICLF"

      character*1 str

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitSuggestedDir$'
     >                   ,0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitSuggestedDir$'
     >                  ,0.,0)

      if (str == 'x' .or. str == 'X') then 
         ppiclf_nbin_dir(1) = 1
      elseif (str == 'y' .or. str == 'Y') then 
        ppiclf_nbin_dir(2) = 1
      elseif (str == 'z' .or. str == 'Z') then 
        if (ppiclf_ndim .lt. 3)
     >  call ppiclf_exittr('Dim must be 3 to use InitSuggestedDir on z$'
     >                   ,0.,ppiclf_ndim)
        ppiclf_nbin_dir(3) = 1
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetNeighborBin
#include "PPICLF"

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
      subroutine ppiclf_solve_NearestNeighbor(i)
#include "PPICLF"

      integer i

      real ydum(PPICLF_LRS), rpropdum(PPICLF_LRP)
      real A(3),B(3),C(3),D(3),AB(3),AC(3)

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

c     do j=1,ppiclf_nwall

c        rnx = ppiclf_wall(1,j)
c        rny = ppiclf_wall(2,j)
c        rnz = 0.0
c        if (ppiclf_ndim .eq. 3) rnz = ppiclf_wall(3,j)
c        rpx = ppiclf_wall(4,j)
c        rpy = ppiclf_wall(5,j)
c        rpz = 0.0
c        if (ppiclf_ndim .eq. 3) rpz = ppiclf_wall(6,j)

c        rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

c        rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
c    >              +rnz*ppiclf_cp_map(3,i)+rd)
c        rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

c        ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
c        ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
c        if (ppiclf_ndim .eq. 3) 
c    >   ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz

c        j_ii = floor((ydum(1)-ppiclf_binb(1))/ppiclf_d2chk(3))
c        j_jj = floor((ydum(2)-ppiclf_binb(3))/ppiclf_d2chk(3))
c        j_kk = floor((ydum(3)-ppiclf_binb(5))/ppiclf_d2chk(3))

c        if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
c        if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
c        if (ppiclf_ndim .eq. 3) then
c        if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
c        endif

c        jp = 0
c        call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
c    >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
c    >                                 ,ydum
c    >                                 ,rpropdum)

c     enddo

      istride = ppiclf_ndim

      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0
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
         if (rflip .gt. 0.0) then
            rnx = -1.*rnx
            rny = -1.*rny
            rnz = -1.*rnz
         endif


         a_sum = 0.0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            endif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            if (rdist .gt. 2.0*ppiclf_d2chk(3)) goto 1511

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0

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
               tri_area = 0.5*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            endif
            a_sum = a_sum + tri_area
         enddo

c        write(6,*) 'Made it', ppiclf_cp_map(1,i),ppiclf_cp_map(2,i)
c    >                       , ppiclf_cp_map(3,i)
c        write(6,*) 'Made it', ydum(1),ydum(2),ydum(3)
c        write(6,*) 'nodes1', ppiclf_wall_c(1,j),ppiclf_wall_c(2,j)
c    >                      , ppiclf_wall_c(3,j)
c        write(6,*) 'nodes2', ppiclf_wall_c(4,j),ppiclf_wall_c(5,j)
c    >                      , ppiclf_wall_c(6,j)
c        write(6,*) 'nodes3', ppiclf_wall_c(7,j),ppiclf_wall_c(8,j)
c    >                      , ppiclf_wall_c(9,j)
c        write(6,*) 'nodes4', ppiclf_wall_c(10,j),ppiclf_wall_c(11,j)
c    >                      , ppiclf_wall_c(12,j)
c        write(6,*) 'normsa', ppiclf_wall_n(1,j),ppiclf_wall_n(2,j)
c    >                      , ppiclf_wall_n(3,j), ppiclf_wall_n(4,j)
c        write(6,*) 'normsf', rnx,rny,rnz
c        write(6,*) 'dist'  , rdist

c        write(6,*) 'Element', j, a_sum, area
         rthresh = 1.10 ! keep it from slipping through crack on edges
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
#include "PPICLF"

      real xp1(*)
      real xp2(*)
      real xp3(*)
      real A(3),B(3),C(3),D(3),AB(3),AC(3)

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitWall$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.,0)

      ppiclf_nwall = ppiclf_nwall + 1 

      if (ppiclf_nwall .gt. PPICLF_LWALL)
     >call ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.,ppiclf_nwall)

      istride = ppiclf_ndim
      a_sum = 0.0
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

         A(1) = (xp1(1) + xp2(1) + xp3(1))/3.0
         A(2) = (xp1(2) + xp2(2) + xp3(2))/3.0
         A(3) = (xp1(3) + xp2(3) + xp3(3))/3.0
      elseif (ppiclf_ndim .eq. 2) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(2)

         A(1) = (xp1(1) + xp2(1))/2.0
         A(2) = (xp1(2) + xp2(2))/2.0
         A(3) = 0.0
      endif

      ! compoute area:
      do k=1,kmax 
         kp = k+1
         if (kp .gt. kmax) kp = kp-kmax ! cycle
         
         kk   = istride*(k-1)
         kkp  = istride*(kp-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = 0.0
         rpx2 = ppiclf_wall_c(kkp+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kkp+2,ppiclf_nwall)
         rpz2 = 0.0

         B(1) = rpx1
         B(2) = rpy1
         B(3) = 0.0
        
         C(1) = rpx2
         C(2) = rpy2
         C(3) = 0.
        
         AB(1) = B(1) - A(1)
         AB(2) = B(2) - A(2)
         AB(3) = 0.0
        
         AC(1) = C(1) - A(1)
         AC(2) = C(2) - A(2)
         AC(3) = 0.0

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
             tri_area = 0.5*AB_MAG*AC_MAG*sin(theta)
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

c     if (ppiclf_nid .eq. 0) then

c        write(6,*) ppiclf_wall_n(1,ppiclf_nwall)
c    >             ,ppiclf_wall_n(2,ppiclf_nwall)
c    >             ,ppiclf_wall_n(3,ppiclf_nwall)
c    >             ,ppiclf_wall_n(4,ppiclf_nwall)
c     endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
#include "PPICLF"

      real xl
      real xr

      if (xl .ge. xr)
     >call ppiclf_exittr('PeriodicX must have xl < xr$',xl,0)

      ppiclf_iperiodic(1) = 0

      ppiclf_xdrange(1,1) = xl
      ppiclf_xdrange(2,1) = xr

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
#include "PPICLF"

      real yl
      real yr

      if (yl .ge. yr)
     >call ppiclf_exittr('PeriodicY must have yl < yr$',yl,0)

      ppiclf_iperiodic(2) = 0

      ppiclf_xdrange(1,2) = yl
      ppiclf_xdrange(2,2) = yr

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
#include "PPICLF"

      real zl
      real zr

      if (zl .ge. zr)
     >call ppiclf_exittr('PeriodicZ must have zl < zr$',zl,0)
      if (ppiclf_ndim .lt. 3)
     >call ppiclf_exittr('Cannot do PeriodicZ if not 3D$',zl,0)

      ppiclf_iperiodic(3) = 0

      ppiclf_xdrange(1,3) = zl
      ppiclf_xdrange(2,3) = zr

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,ngrid)
#include "PPICLF"

      real    filt
      real    alpha
      integer ngrid

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.,0)

      ppiclf_filter = filt
      ppiclf_ngrids = ngrid
      ppiclf_alpha  = alpha 

      rsig             = ppiclf_filter/(2.*sqrt(2.*log(2.)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_alpha))

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT      = .true.
      PPICLF_LFILTGAUSS = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitBoxFilter(filt,ngrid)
#include "PPICLF"

      real    filt
      integer ngrid

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.,0)

      ppiclf_filter = filt
      ppiclf_ngrids = ngrid

      ppiclf_d2chk(2)  = filt/2.0

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT    = .true.
      PPICLF_LFILTBOX = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetParticleTag
#include "PPICLF"

      do i=1,ppiclf_npart
         if (ppiclf_iprop(5,i) .eq. -1) ppiclf_iprop(5,i) = ppiclf_nid 
         if (ppiclf_iprop(6,i) .eq. -1) ppiclf_iprop(6,i) = ppiclf_cycle
         if (ppiclf_iprop(7,i) .eq. -1) ppiclf_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParticleTag
#include "PPICLF"

      if (.not. PPICLF_RESTART) then
         do i=1,PPICLF_LPART
            ppiclf_iprop(5,i) = -1
            ppiclf_iprop(6,i) = -1
            ppiclf_iprop(7,i) = -1
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time
     >                                         ,y,ydot)
#include "PPICLF"

      integer istep
      integer iostep
      real    dt
      real    time
      real    y(*)
      real    ydot(*)

      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      if (ppiclf_imethod .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(time,y,ydot)
      if (ppiclf_imethod .eq. 2) 
     >   call ppiclf_solve_IntegrateFwdEuler(time,y,ydot)

      ! output files
      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0) then

         ! already wrote initial conditions
         if (ppiclf_cycle .ne. 0) then
            call ppiclf_io_WriteParticleVTU('')
            call ppiclf_io_WriteBinVTU('')
         endif

         if (ppiclf_lsubbin)
     >      call ppiclf_io_WriteSubBinVTU('')

         call ppiclf_io_OutputDiagAll
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3(time_,y,ydot)
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      ndum = PPICLF_NPART*PPICLF_LRS

      ! save stage 1 solution
      do i=1,ndum
         ppiclf_y1(i) = y(i)
      enddo

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_solve_SetYdot(time_,y,ydot)

         ! rk3 integrate
         do i=1,ndum
            y(i) =  ppiclf_rk3coef(1,istage)*ppiclf_y1 (i)
     >            + ppiclf_rk3coef(2,istage)*y      (i)
     >            + ppiclf_rk3coef(3,istage)*ydot   (i)
         enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateFwdEuler(time_,y,ydot)
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      ndum = PPICLF_NPART*PPICLF_LRS

      ! save stage 1 solution
      do i=1,ndum
         ppiclf_y1(i) = y(i)
      enddo

      ! evaluate ydot
      call ppiclf_solve_SetYdot(time_,y,ydot)

      ! rk3 integrate
      do i=1,ndum
         y(i) = y(i) + ppiclf_dt*ydot(i) 
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetYdot(time_,y,ydot)
#include "PPICLF"

      real time_
      real y(*)
      real ydot(*)

      call ppiclf_solve_InitSolve
      call ppiclf_user_SetYdot(time_,y,ydot)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSolve
#include "PPICLF"

      call ppiclf_solve_RemoveParticle
      call ppiclf_comm_CreateBin
      call ppiclf_comm_FindParticle
      call ppiclf_comm_MoveParticle
      if (ppiclf_overlap) call ppiclf_comm_MapOverlapMesh

      if (ppiclf_lsubsubbin .or. ppiclf_lproj) then
         call ppiclf_comm_CreateGhost
         call ppiclf_comm_MoveGhost
      endif

      if (ppiclf_lproj .and. ppiclf_overlap) 
     >   call ppiclf_solve_ProjectParticleGrid
      if (ppiclf_lsubsubbin) 
     >   call ppiclf_solve_SetNeighborBin

      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ppiclf_ydotc(j,i) = 0.0
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitInterp
#include "PPICLF"

      if (.not.ppiclf_overlap)
     >call ppiclf_exittr('Cannot interpolate unless overlap grid$',0.,0)
      if (PPICLF_INTERP .eq. 0)
     >call ppiclf_exittr(
     >     'No specified interpolated fields, set PPICLF_LRP_INT$',0.,0)

      PPICLF_INT_ICNT = 0

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
      call fgslib_findpts_setup(ppiclf_fp_hndl
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
      subroutine ppiclf_solve_InterpField(jp,infld)
#include "PPICLF"

      real infld(*)

      if (ppiclf_int_icnt .eq. -1)
     >call ppiclf_exittr('Call InitInterp before InterpField$',0.,0)

      PPICLF_INT_ICNT = PPICLF_INT_ICNT + 1

      if (PPICLF_INT_ICNT .gt. PPICLF_LRP_INT)
     >   call ppiclf_exittr('Interpolating too many fields$'
     >                     ,0.0,PPICLF_INT_ICNT)
      if (jp .le. 0 .or. jp .gt. PPICLF_LRP)
     >   call ppiclf_exittr('Invalid particle array interp. location$'
     >                     ,0.0,jp)

      PPICLF_INT_MAP(PPICLF_INT_ICNT) = jp

      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltbbb
         iee = ppiclf_er_mapc(1,ie)
         j = (iee-1)*n + 1
         call ppiclf_copy(ppiclf_int_fld(1,1,1,PPICLF_INT_ICNT,ie)
     >                   ,infld(j),n)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_FinalizeInterp
#include "PPICLF"

      REAL FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)

      integer nkey(2)

      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltbbb
     >      ,PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld
     >      ,nrr,njj)
      call fgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltbbb
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld,nrr,nkey,2)

      ! find which cell particle is in locally
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim .eq. 3)
     >iz = 3

      call fgslib_findpts(PPICLF_FP_HNDL           !   call fgslib_findpts( ihndl,
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
         call fgslib_findpts_eval_local( PPICLF_FP_HNDL
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
      call fgslib_findpts_free(PPICLF_FP_HNDL)

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff(dt)
#include "PPICLF"

      real dt

      ppiclf_rk3coef(1,1) = 0.0
      ppiclf_rk3coef(2,1) = 1.0 
      ppiclf_rk3coef(3,1) = dt
      ppiclf_rk3coef(1,2) = 3.0/4.0
      ppiclf_rk3coef(2,2) = 1.0/4.0 
      ppiclf_rk3coef(3,2) = dt/4.0 
      ppiclf_rk3coef(1,3) = 1.0/3.0
      ppiclf_rk3coef(2,3) = 2.0/3.0 
      ppiclf_rk3coef(3,3) = dt*2.0/3.0 

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_RemoveParticle
#include "PPICLF"

      integer in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >                                   iperiodicz,ndim

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
         do j=0,ndim-1
            jchk = jj(j+1)
            if (ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
                   ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
                  goto 1512
                endif
            endif
            if (ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1))then
               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
                   ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) +
     &                     abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
                   ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
     &                     abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (ppiclf_iprop(1,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1512 continue
         enddo
      enddo

      ic = 0
      do i=1,ppiclf_npart
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
      ppiclf_npart = ic

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleGrid
#include "PPICLF"

      real    multfci

      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      integer ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp

      logical partl, if3d

      integer nkey(2)

      real pi

      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      nxyzdum = nxyz*PPICLF_LRP_PRO*PPICLF_LEE
      do i=1,nxyzdum
         ppiclf_pro_fldb(i,1,1,1,1) = 0.0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d) ppiclf_jzgp  = 3

      rdum = 0.0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d+0)
         rdum   = 1./(-2.*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         multfci = 1.0/ppiclf_filter**2
         if (if3d) multfci = multfci/ppiclf_filter
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
                    
         iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
         iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
         if (if3d)
     >   iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
         iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp

      do ip=1,ndum
         iip      = iproj(1,ip)
         jjp      = iproj(2,ip)
         if (if3d)
     >   kkp      = iproj(3,ip)
         ndumdum  = iproj(4,ip)

         ilow  = iip-1
         ihigh = iip+1
         jlow  = jjp-1
         jhigh = jjp+1
         if (if3d) then
            klow  = kkp-1
            khigh = kkp+1
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
     >                (ppiclf_xm1b(i,j,k,2,ie) - rproj(3,ip))**2
            if(if3d) rdist2 = rdist2 +
     >                (ppiclf_xm1b(i,j,k,3,ie) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = 1.0
            if (ppiclf_lfiltgauss)
     >         rexp = exp(rdist2*rproj(1,ip))
            
            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                         ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                       + rproj(j1,ip)*rexp
            enddo
         enddo
         enddo
         enddo
         enddo
      enddo

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
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc,
     >   PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)
      call fgslib_crystal_tuple_sort    (ppiclf_cr_hndl,neltbc
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,nkey,2)

      ! add the fields from the bins to ptw array
      nlxyzep = nxyz*PPICLF_LEE*PPICLF_LRP_PRO
      do i=1,nlxyzep
         PPICLF_PRO_FLD(i,1,1,1,1) = 0.0
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
#include "PPICLF"

      real    multfci
      real    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer iproj(4,PPICLF_LPART+PPICLF_LPART_GP)

      integer ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp

      logical if3d

      real pi

      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d)
     >ppiclf_jzgp  = 3

      rdum = 0.0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.*sqrt(2.*log(2.)))
         multfci = 1./(sqrt(2.*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d+0)
         rdum   = 1./(-2.*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         multfci = 1.0/ppiclf_filter**2
         if (if3d) multfci = multfci/ppiclf_filter
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
         idum = floor(ppiclf_filter/2.0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1
         jdum = floor(ppiclf_filter/2.0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1
         if (if3d)
     >   kdum = floor(ppiclf_filter/2.0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_alpha)/log(2.0)))+1
      endif

      if (ppiclf_lfiltbox) then
         idum = ngrids/2+1
         jdum = ngrids/2+1
         if (if3d)
     >   kdum = ngrids/2+1
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

         do k=kl,kr
         do j=jl,jr
         do i=il,ir
            rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
            if(if3d) rdist2 = rdist2 +
     >                (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = 1.0
            if (ppiclf_lfiltgauss)
     >         rexp = exp(rdist2*rproj(1,ip))

            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_grid_fld(i,j,k,jj) = 
     >                         ppiclf_grid_fld(i,j,k,jj) 
     >                       + rproj(j1,ip)*rexp
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
