!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitMPI(comm,id,np)
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Input: 
!
      integer*4 comm
      integer*4 id
      integer*4 np
!
      if (PPICLF_LINIT .or. PPICLF_LFILT .or. PPICLF_OVERLAP)
     >   call ppiclf_exittr('InitMPI must be called first$',0.0,0)

      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      call ppiclf_prints('   *Begin InitCrystal$')
         call ppiclf_comm_InitCrystal
      call ppiclf_prints('    End InitCrystal$')

      PPICLF_LCOMM = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitFindptsDum
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      real*8 tol, bb_t, npt_max
      integer*4 n, ie
      real*4 xm1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE), 
     >       ym1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE), 
     >       zm1(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
      common /ppiclf_tmp_grid/ xm1, ym1, zm1
!
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
         call ppiclf_copy(xm1(1,1,1,ie),ppiclf_xm1b(1,1,1,1,ie),n)
         call ppiclf_copy(ym1(1,1,1,ie),ppiclf_xm1b(1,1,1,2,ie),n)
         call ppiclf_copy(zm1(1,1,1,ie),ppiclf_xm1b(1,1,1,3,ie),n)
      enddo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128

      call fgslib_findpts_setup(ppiclf_fp_hndl
     >                         ,ppiclf_comm_nid
     >                         ,1 ! only 1 rank on this comm
     >                         ,ppiclf_ndim
     >                         ,xm1
     >                         ,ym1
     >                         ,zm1
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_neltb
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ppiclf_neltb+2
     >                         ,ppiclf_neltb+2
     >                         ,npt_max
     >                         ,tol)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitCrystal
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
      call fgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateBin
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer*4  ifac(3), icount(3)
      real*8     d2new(3)
      integer*4 ix, iy, iz, iperiodicx, iperiodicy, iperiodicz, ndim,
     >          npt_total, j, i, nmax, nbb, idum, jdum, kdum, nbin
      real*8 xmin, ymin, zmin, xmax, ymax, zmax, rduml, rdumr
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      real*8 ppiclf_glmin,ppiclf_glmax
      external ppiclf_glmin,ppiclf_glmax
!

! face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq. 3)
     >iz = 3

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
      ndim       = ppiclf_ndim
         
      ppiclf_d2chk(1) = max(ppiclf_d2chk(2),ppiclf_d2chk(3))

      ! binning requires > 1 global particle. This takes care of 
      ! single particle case
      npt_total = ppiclf_iglsum(ppiclf_npart,1)
      if (npt_total .eq. 1) then
      if (.not. ppiclf_lproj .and. .not. ppiclf_lsubsubbin) 
     >ppiclf_d2chk(1) = 1E-3
      endif

      ! compute binb
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      do i=1,ppiclf_npart
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(1)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(1)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         if (ppiclf_ndim .eq. 3) then
            rduml = ppiclf_y(iz,i) - ppiclf_d2chk(1)
            rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(1)
            if (rduml .lt. zmin) zmin = rduml
            if (rdumr .gt. zmax) zmax = rdumr
         endif
      enddo

      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = 0.0
      ppiclf_binb(6) = 0.0
      if(ppiclf_ndim .gt. 2) ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      if(ppiclf_ndim .gt. 2) ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      if (ppiclf_xdrange(2,1) .lt. ppiclf_binb(2) .or.
     >    ppiclf_xdrange(1,1) .gt. ppiclf_binb(1)) then
         ppiclf_binb(1) = ppiclf_xdrange(1,1)
         ppiclf_binb(2) = ppiclf_xdrange(2,1)
      endif

      if (ppiclf_xdrange(2,2) .lt. ppiclf_binb(4) .or.
     >    ppiclf_xdrange(1,2) .gt. ppiclf_binb(3)) then
         ppiclf_binb(3) = ppiclf_xdrange(1,2)
         ppiclf_binb(4) = ppiclf_xdrange(2,2)
      endif

      if (ppiclf_ndim .gt. 2) then
      if (ppiclf_xdrange(2,3) .lt. ppiclf_binb(6) .or.
     >    ppiclf_xdrange(1,3) .gt. ppiclf_binb(5)) then
         ppiclf_binb(5) = ppiclf_xdrange(1,3)
         ppiclf_binb(6) = ppiclf_xdrange(2,3)
      endif
      endif

      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      icount(1) = 0
      icount(2) = 0
      icount(3) = 0
      d2new(1) = ppiclf_d2chk(1)
      d2new(2) = ppiclf_d2chk(1)
      d2new(3) = ppiclf_d2chk(1)

      ! first check if suggested direction 
      j = -1
      nmax = int(1E8)
      if (ppiclf_nbin_dir(1) .eq. 1) j = 0
      if (ppiclf_nbin_dir(2) .eq. 1) j = 1
      if (ppiclf_nbin_dir(3) .eq. 1) j = 2

      if (j .ge. 0) then
         do i=1,nmax
            ifac(j+1) = ifac(j+1) + 1
            d2new(j+1) = (ppiclf_binb(2+2*j) - ppiclf_binb(1+2*j))/
     >                    ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)
           
            if( nbb .gt. ppiclf_np ) then
               icount(j+1) = 1
               ifac(j+1) = ifac(j+1) - 1
               d2new(j+1) = (ppiclf_binb(2+2*j) -ppiclf_binb(1+2*j))/
     >                       ifac(j+1)
            endif
            if (icount(j+1) .gt. 0) then
               exit
            endif
         enddo
      endif


      ! then do every direction
      do i=1,nmax
         do j=0,ndim-1
            ifac(j+1) = ifac(j+1) + 1
            d2new(j+1) = (ppiclf_binb(2+2*j) - ppiclf_binb(1+2*j))/
     >                    ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)
         
            if( nbb .gt. ppiclf_np ) then
c              icount(j+1) = icount(j+1) + 1
               icount(j+1) = 1
               ifac(j+1) = ifac(j+1) - 1
               d2new(j+1) = (ppiclf_binb(2+2*j) -ppiclf_binb(1+2*j))/
     >                       ifac(j+1)
            endif
         enddo
         if (icount(1) .gt. 0) then
         if (icount(2) .gt. 0) then
         if (icount(3) .gt. 0 .or. ndim .lt. 3) then
            exit
         endif
         endif
         endif
      enddo

! -------------------------------------------------------
! SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! how many spacings in each direction
      ppiclf_ndxgp = floor( (ppiclf_binb(2) - ppiclf_binb(1))/d2new(1))
      ppiclf_ndygp = floor( (ppiclf_binb(4) - ppiclf_binb(3))/d2new(2))
      ppiclf_ndzgp = 1
      if (ppiclf_ndim .gt. 2) ppiclf_ndzgp = 
     >                floor( (ppiclf_binb(6) - ppiclf_binb(5))/d2new(3))
      
      ! grid spacing for that many spacings
      ppiclf_rdxgp = (ppiclf_binb(2) -ppiclf_binb(1))/real(ppiclf_ndxgp)
      ppiclf_rdygp = (ppiclf_binb(4) -ppiclf_binb(3))/real(ppiclf_ndygp)
      ppiclf_rdzgp = 1.
      if (ppiclf_ndim .gt. 2) 
     >ppiclf_rdzgp = (ppiclf_binb(6) -ppiclf_binb(5))/real(ppiclf_ndzgp)

      nbin = ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp

!     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_ndxgp)
         jdum = modulo(ppiclf_nid/ppiclf_ndxgp,ppiclf_ndygp)
         kdum = ppiclf_nid/(ppiclf_ndxgp*ppiclf_ndygp)
         if (ppiclf_ndim .lt. 3) kdum = 0
         ppiclf_binx(1,1) = ppiclf_binb(1) + idum    *ppiclf_rdxgp
         ppiclf_binx(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_rdxgp
         ppiclf_biny(1,1) = ppiclf_binb(3) + jdum    *ppiclf_rdygp
         ppiclf_biny(2,1) = ppiclf_binb(3) + (jdum+1)*ppiclf_rdygp
         ppiclf_binz(1,1) = 0.0
         ppiclf_binz(2,1) = 0.0
         if (ppiclf_ndim .gt. 2) then
            ppiclf_binz(1,1) = ppiclf_binb(5) + kdum    *ppiclf_rdzgp
            ppiclf_binz(2,1) = ppiclf_binb(5) + (kdum+1)*ppiclf_rdzgp
         endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateSubBin
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      integer*4 nbin, idum, jdum, kdum, ndumx, ndumy, itmp, jtmp, ktmp,
     >          i, j, k
!

      nbin = ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp

c     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_ndxgp)
         jdum = modulo(ppiclf_nid/ppiclf_ndxgp,ppiclf_ndygp)
         kdum = ppiclf_nid/(ppiclf_ndxgp*ppiclf_ndygp)
         if (ppiclf_ndim .lt. 3) kdum = 0
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_rdxgp/ppiclf_filter) + 1 + 1
         ppiclf_by = floor(ppiclf_rdygp/ppiclf_filter) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_rdzgp/ppiclf_filter) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_ngrids
         ppiclf_by = ppiclf_by*ppiclf_ngrids
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_ngrids

         if (ppiclf_bx .gt. PPICLF_BX1)
     >      call ppiclf_exittr('Increase PPICLF_BX1$',0.,ppiclf_bx)
         if (ppiclf_by .gt. PPICLF_BY1)
     >      call ppiclf_exittr('Increase PPICLF_BY1$',0.,ppiclf_by)
         if (ppiclf_bz .gt. PPICLF_BZ1)
     >      call ppiclf_exittr('Increase PPICLF_BZ1$',0.,ppiclf_bz)

         ppiclf_rdx = ppiclf_rdxgp/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_rdygp/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_rdz = ppiclf_rdzgp/(ppiclf_bz-1)

         ndumx = ppiclf_ndxgp*(ppiclf_bx-1) + 1
         ndumy = ppiclf_ndygp*(ppiclf_by-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = (ppiclf_binx(1,1) +
     >                                  (i-1)*ppiclf_rdx)
            ppiclf_grid_y(i,j,k) = (ppiclf_biny(1,1) +
     >                                  (j-1)*ppiclf_rdy)
            ppiclf_grid_z(i,j,k) = (ppiclf_binz(1,1) +
     >                                  (k-1)*ppiclf_rdz)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         enddo
         enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MapOverlapMesh
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
      include 'mpif.h'
!
! Internal:
!
      integer*4 icalld
      save      icalld
      data      icalld /0/
      integer*4 nkey(2), i, j, k, ie, iee, ii, jj, kk, ndum, nrank,
     >          nl, nii, njj, nrr, ilow, jlow, klow, nxyz, il,
     >          ihigh, jhigh, khigh, ierr
      real*8 rxval, ryval, rzval
      logical partl
      real*8 ppiclf_vlmin, ppiclf_vlmax
      external ppiclf_vlmin, ppiclf_vlmax
!

      ! see which bins are in which elements
      ppiclf_neltb = 0
      do ie=1,ppiclf_nee
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         rxval = ppiclf_xm1bs(i,j,k,1,ie)
         ryval = ppiclf_xm1bs(i,j,k,2,ie)
         rzval = 0.
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1bs(i,j,k,3,ie)

         if (rxval .gt. ppiclf_binb(2)) goto 1233
         if (rxval .lt. ppiclf_binb(1)) goto 1233
         if (ryval .gt. ppiclf_binb(4)) goto 1233
         if (ryval .lt. ppiclf_binb(3)) goto 1233
         if (ppiclf_ndim.gt.2 .and. rzval .gt. ppiclf_binb(6)) 
     >      goto 1233
         if (ppiclf_ndim.gt.2 .and. rzval .lt. ppiclf_binb(5))
     >      goto 1233

         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_rdxgp) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_rdygp) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_rdzgp) 
         if (ppiclf_ndim.lt.3) kk = 0
          if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
          if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
          if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
          if (ii .eq. -1) ii = 0
          if (jj .eq. -1) jj = 0
          if (kk .eq. -1) kk = 0
         ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk
         nrank = ndum

c        if (ii .eq. ppiclf_ndxgp-1) write(6,*) rxval,ppiclf_ndxgp

c        if (ppiclf_binb(3) .gt. ryval)
c    >      write(6,*) ryval,ppiclf_binb(3),jj
         if (ii .lt. 0 .or. ii .gt. ppiclf_ndxgp-1) then
c           write(6,*) 'Bounds here:',ppiclf_binb
c           write(6,*) 'Failed here:',rxval,ryval,rzval
            goto 1233
         endif
         if (jj .lt. 0 .or. jj .gt. ppiclf_ndygp-1) then
c           write(6,*) 'Bounds here:',ppiclf_binb
c           write(6,*) 'Failed here:',rxval,ryval,rzval
            goto 1233
         endif
         if (kk .lt. 0 .or. kk .gt. ppiclf_ndzgp-1) then
c           write(6,*) 'Bounds here:',ppiclf_binb
c           write(6,*) 'Failed here:',rxval,ryval,rzval
            goto 1233
         endif

         ppiclf_neltb = ppiclf_neltb + 1
         if(ppiclf_neltb .gt. PPICLF_LEE) then
           call ppiclf_exittr('Increase PPICLF_LEE$',0.,ppiclf_neltb)
         endif

         ppiclf_er_map(1,ppiclf_neltb) = ie
         ppiclf_er_map(2,ppiclf_neltb) = ppiclf_nid
         ppiclf_er_map(3,ppiclf_neltb) = ndum
         ppiclf_er_map(4,ppiclf_neltb) = nrank
         ppiclf_er_map(5,ppiclf_neltb) = nrank
         ppiclf_er_map(6,ppiclf_neltb) = nrank

         if (ppiclf_neltb .gt. 1) then
         do il=1,ppiclf_neltb-1
            if (ppiclf_er_map(1,il) .eq. ie) then
            if (ppiclf_er_map(4,il) .eq. nrank) then
               ppiclf_neltb = ppiclf_neltb - 1
               goto 1233
            endif
            endif
         enddo
         endif
 1233 continue
      enddo
      enddo
      enddo
      enddo

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
       iee = ppiclf_er_map(1,ie)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,1,ie)
     >                 ,ppiclf_xm1bs(1,1,1,1,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,2,ie)
     >                 ,ppiclf_xm1bs(1,1,1,2,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,3,ie)
     >                 ,ppiclf_xm1bs(1,1,1,3,iee),nxyz)
      enddo

      ppiclf_neltbb = ppiclf_neltb
      do ie=1,ppiclf_neltbb
         call ppiclf_icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie)
     >             ,PPICLF_LRMAX)
      enddo


      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*3
      nkey(1) = 2
      nkey(2) = 1
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltb
     >       ,PPICLF_LEE,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
      call fgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltb
     >       ,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,2)


      do ie=1,ppiclf_neltb
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         rxval = ppiclf_xm1b(i,j,k,1,ie)
         ryval = ppiclf_xm1b(i,j,k,2,ie)
         rzval = 0.
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1b(i,j,k,3,ie)
         
         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_rdxgp) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_rdygp) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_rdzgp) 
         if (ppiclf_ndim.eq.2) kk = 0
          if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
          if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
          if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
          if (ii .eq. -1) ii = 0
          if (jj .eq. -1) jj = 0
          if (kk .eq. -1) kk = 0
         ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk

         ppiclf_modgp(i,j,k,ie,1) = ii
         ppiclf_modgp(i,j,k,ie,2) = jj
         ppiclf_modgp(i,j,k,ie,3) = kk
         ppiclf_modgp(i,j,k,ie,4) = ndum
   
      enddo
      enddo
      enddo
      enddo

      do ie=1,ppiclf_neltb
         ppiclf_xerange(1,1,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(2,1,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(1,2,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(2,2,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(1,3,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         ppiclf_xerange(2,3,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)

         ilow  = 
     >     floor((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/ppiclf_rdxgp)
         ihigh = 
     >     floor((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/ppiclf_rdxgp)
         jlow  = 
     >     floor((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/ppiclf_rdygp)
         jhigh = 
     >     floor((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/ppiclf_rdygp)
         klow  = 
     >     floor((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/ppiclf_rdzgp)
         khigh = 
     >     floor((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/ppiclf_rdzgp)
         if (ppiclf_ndim.lt.3) then
            klow = 0
            khigh = 0
         endif

         ppiclf_el_map(1,ie) = ilow  + ppiclf_ndxgp*jlow  
     >                            + ppiclf_ndxgp*ppiclf_ndygp*klow
         ppiclf_el_map(2,ie) = ihigh + ppiclf_ndxgp*jhigh 
     >                            + ppiclf_ndxgp*ppiclf_ndygp*khigh
         ppiclf_el_map(3,ie) = ilow
         ppiclf_el_map(4,ie) = ihigh
         ppiclf_el_map(5,ie) = jlow
         ppiclf_el_map(6,ie) = jhigh
         ppiclf_el_map(7,ie) = klow
         ppiclf_el_map(8,ie) = khigh
      enddo

      if (icalld .eq. 0) then 

         icalld = icalld + 1

         call ppiclf_prints('   *Begin mpi_comm_split$')
            call mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr)
         call ppiclf_prints('    End mpi_comm_split$')

         call ppiclf_io_OutputDiagGrid
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Input:
!
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8    xgrid(*)
      real*8    ygrid(*)
      real*8    zgrid(*)
!
! External:
!
      integer*4 nxyz, i, j, ie
!

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitOverlap$',0.,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.,0)

      if (ncell .gt. PPICLF_LEE .or. ncell .lt. 0) 
     >   call ppiclf_exittr('Increase LEE in InitOverlap$',0.,ncell)
      if (lx1 .ne. PPICLF_LEX) 
     >   call ppiclf_exittr('LX1 != LEX in InitOverlap$',0.,ncell)
      if (ly1 .ne. PPICLF_LEY)
     >   call ppiclf_exittr('LY1 != LEY in InitOverlap$',0.,ncell)
      if (lz1 .ne. PPICLF_LEZ)
     >   call ppiclf_exittr('LZ1 != LEZ in InitOverlap$',0.,ncell)

      ppiclf_nee = ncell
      nxyz       = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      do ie=1,ppiclf_nee
      do i=1,nxyz
         j = (ie-1)*nxyz + i
         ppiclf_xm1bs(i,1,1,1,ie) = xgrid(j)
         ppiclf_xm1bs(i,1,1,2,ie) = ygrid(j)
         ppiclf_xm1bs(i,1,1,3,ie) = zgrid(j)
      enddo
      enddo

      call ppiclf_comm_MapOverlapMesh

      ppiclf_overlap = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_FindParticle
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      integer*4 ix, iy, iz, i, ii, jj, kk, ndum, nrank
!
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq.3)
     >iz = 3

      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ii    = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_rdxgp) 
         jj    = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_rdygp) 
         kk    = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_rdzgp) 
         if (ppiclf_ndim .lt. 3) kk = 0
c        if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
c        if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
c        if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
         ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk
         nrank = ndum

         ppiclf_iprop(8,i)  = ii
         ppiclf_iprop(9,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(3,i)  = nrank ! where particle is actually moved
         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveParticle
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      logical partl    
      integer*4 lrf
      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2)
      real*8 rwork(lrf,PPICLF_LPART)
      integer*4 i, ic, j0
!

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
      enddo

      j0 = 4
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart,PPICLF_LPART
     >                                  ,ppiclf_iprop,PPICLF_LIP
     >                                  ,partl,0
     >                                  ,rwork,lrf
     >                                  ,j0)

      if (ppiclf_npart .gt. PPICLF_LPART .or. ppiclf_npart .lt. 0)
     >   call ppiclf_exittr('Increase LPART$',0.0,ppiclf_npart)

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
      enddo
        
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateGhost
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         rxval = ppiclf_cp_map(1,ip)
         ryval = ppiclf_cp_map(2,ip)
         rzval = 0.
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip)

         iip    = ppiclf_iprop(8,ip)
         jjp    = ppiclf_iprop(9,ip)
         kkp    = ppiclf_iprop(10,ip)

         rxl = ppiclf_binb(1) + ppiclf_rdxgp*iip
         rxr = rxl + ppiclf_rdxgp
         ryl = ppiclf_binb(3) + ppiclf_rdygp*jjp
         ryr = ryl + ppiclf_rdygp
         rzl = 0.0
         rzr = 0.0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_rdzgp*kkp
            rzr = rzl + ppiclf_rdzgp
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_ndxgp*jjg 
     >                  + ppiclf_ndxgp*ppiclf_ndygp*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_ndxgp*jjg 
     >                  + ppiclf_ndxgp*ppiclf_ndygp*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_ndxgp*jjg 
     >                  + ppiclf_ndxgp*ppiclf_ndygp*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Input:
!
      real*8 rxdrng(3)
      integer*4 iadd(3)
!
! Input/Output:
!
      real*8 rxnew(3)
!
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .ge. ppiclf_ndxgp) then
         rxnew(1) = rxnew(1) - rxdrng(1)
         goto 123
      endif
      endif
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .lt. 0) then
         rxnew(1) = rxnew(1) + rxdrng(1)
         goto 123
      endif
      endif

  123 continue    
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .ge. ppiclf_ndygp) then
         rxnew(2) = rxnew(2) - rxdrng(2)
         goto 124
      endif
      endif
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .lt. 0) then
         rxnew(2) = rxnew(2) + rxdrng(2)
         goto 124
      endif
      endif
  124 continue

      if (ppiclf_ndim .gt. 2) then
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .ge. ppiclf_ndzgp) then
            rxnew(3) = rxnew(3) - rxdrng(3)
            goto 125
         endif
         endif
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .lt. 0) then
            rxnew(3) = rxnew(3) + rxdrng(3)
            goto 125
         endif
         endif
      endif
  125 continue

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveGhost
!
      implicit none
!
#include "PPICLF.h"
#include "PPICLF"      
!
! Internal:
!
      logical partl         
!
      call fgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
