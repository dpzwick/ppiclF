!-----------------------------------------------------------------------
      subroutine ppiclf_comm_setup
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

c     call ppiclf_comm_interp_setup(i_fp_hndl,0.0,idum,ppiclf_nelt)
      call fgslib_crystal_setup(i_cr_hndl,ppiclf_comm,ppiclf_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_interp_setup(ih,tolin,nmsh,nelm)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      ! need to add xm1, etc in PPICLF for this case

c     common /intp_h/ ih_intp(2,1)
c     common /intp/   tol

c     data ihcounter /0/
c     save ihcounter

c     real xmi, ymi, zmi
c     common /CBXMI/ xmi(PPICLF_LX1*PPICLF_LY1*PPICLF_LZ1*PPICLF_LELT),
c    $               ymi(PPICLF_LX1*PPICLF_LY1*PPICLF_LZ1*PPICLF_LELT),
c    $               zmi(PPICLF_LX1*PPICLF_LY1*PPICLF_LZ1*PPICLF_LELT)
c
c     real w(2*PPICLF_LX1*PPICLF_LY1*PPICLF_LZ1)
c     logical if3d

c     tol = max(5e-13,tolin)
c     npt_max = 128
c     bb_t    = 0.01

c     ! fix below
c     if3d = .true.

c     if (nmsh.ge.1 .and. nmsh.lt.PPICLF_LX1-1) then
c        nxi = nmsh+1
c        nyi = nxi
c        nzi = 1
c        if (if3d) nzi = nxi
c        do ie = 1,nelm
c          call map_m_to_n(xmi((ie-1)*nxi*nyi*nzi + 1),nxi,xm1(1,1,1,ie)
c    $                     ,lx1,if3d,w,size(w))
c          call map_m_to_n(ymi((ie-1)*nxi*nyi*nzi + 1),nyi,ym1(1,1,1,ie)
c    $                     ,ly1,if3d,w,size(w))
c          if (if3d) 
c    $     call map_m_to_n(zmi((ie-1)*nxi*nyi*nzi + 1),nzi,zm1(1,1,1,ie)
c    $                     ,lz1,if3d,w,size(w))
c        enddo

c        n = nelm*nxi*nyi*nzi 
c        call fgslib_findpts_setup(ih_intp1,nekcomm,npp,ldim,
c    &                             xm1,ym1,zm1,nx1,ny1,nz1,
c    &                             nelm,nx1,ny1,nz1,bb_t,nelm+2,nelm+2,
c    &                             npt_max,tol)

c        call fgslib_findpts_setup(ih_intp2,nekcomm,npp,ldim,
c    $                             xmi,ymi,zmi,nxi,nyi,nzi,
c    $                             nelm,2*nxi,2*nyi,2*nzi,bb_t,n,n,
c    $                             npt_max,tol)
c     else
c        n = nelm*nx1*ny1*nz1
c        call fgslib_findpts_setup(ih_intp1,nekcomm,npp,ldim,
c    &                             xm1,ym1,zm1,nx1,ny1,nz1,
c    &                             nelm,2*nx1,2*ny1,2*nz1,bb_t,n,n,
c    &                             npt_max,tol)
c        ih_intp2 = ih_intp1
c     endif

c     ihcounter = ihcounter + 1 
c     ih = ihcounter
c     if (ih .gt. INTP_HMAX)
c    $   call exitti('Maximum number of handles exceeded!$',INTP_HMAX)
c     ih_intp(1,ih) = ih_intp1
c     ih_intp(2,ih) = ih_intp2

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_bin_setup
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      integer  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer  ifac(3), icount(3)
      real     d2new(3)
      logical  partl

      real                      ppiclf_xerange(2,3,ppiclf_lbmax)
      common /ppiclf_elementrange/ ppiclf_xerange

c     face, edge, and corner number, x,y,z are all inline, so stride=3
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

      if (ppiclf_rparam(12) .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 3

      iperiodicx = int(ppiclf_rparam(8))
      iperiodicy = int(ppiclf_rparam(9))
      iperiodicz = int(ppiclf_rparam(10))
      ndim       = int(ppiclf_rparam(12))

      ! compute binb
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      do i=1,ppiclf_npart
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(2)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(2)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(2)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(2)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         rduml = ppiclf_y(iz,i) - ppiclf_d2chk(2)
         rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(2)
         if (rduml .lt. zmin) zmin = rduml
         if (rdumr .gt. zmax) zmax = rdumr
      enddo

      ppiclf_binb(1) = glmin(xmin,1)
      ppiclf_binb(2) = glmax(xmax,1)
      ppiclf_binb(3) = glmin(ymin,1)
      ppiclf_binb(4) = glmax(ymax,1)
      ppiclf_binb(5) = 0.0
      ppiclf_binb(6) = 0.0
      if(ppiclf_rparam(12) .gt. 2) ppiclf_binb(5) = glmin(zmin,1)
      if(ppiclf_rparam(12) .gt. 2) ppiclf_binb(6) = glmax(zmax,1)


      ! comment for now.....
c     ppiclf_binb(1) = max(ppiclf_binb(1),ppiclf_xdrange(1,1))
c     ppiclf_binb(2) = min(ppiclf_binb(2),ppiclf_xdrange(2,1))
c     ppiclf_binb(3) = max(ppiclf_binb(3),ppiclf_xdrange(1,2))
c     ppiclf_binb(4) = min(ppiclf_binb(4),ppiclf_xdrange(2,2))
c     if(ppiclf_rparam(12) .gt. 2) ppiclf_binb(5) = 
c    >                           max(ppiclf_binb(5),ppiclf_xdrange(1,3))
c     if(ppiclf_rparam(12) .gt. 2) ppiclf_binb(6) = 
c    >                           min(ppiclf_binb(6),ppiclf_xdrange(2,3))
c     if (iperiodicx .eq. 0) then
c        ppiclf_binb(1) = ppiclf_xdrange(1,1)
c        ppiclf_binb(2) = ppiclf_xdrange(2,1)
c     endif
c     if (iperiodicy .eq. 0) then
c        ppiclf_binb(3) = ppiclf_xdrange(1,2)
c        ppiclf_binb(4) = ppiclf_xdrange(2,2)
c     endif
c     if (iperiodicz .eq. 0 .and. ppiclf_rparam(12) .gt. 2) then
c        ppiclf_binb(5) = ppiclf_xdrange(1,3)
c        ppiclf_binb(6) = ppiclf_xdrange(2,3)
c     endif

      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      icount(1) = 0
      icount(2) = 0
      icount(3) = 0
      d2new(1) = ppiclf_d2chk(2)
      d2new(2) = ppiclf_d2chk(2)
      d2new(3) = ppiclf_d2chk(2)

      ppiclf_ndxgp = floor( (ppiclf_binb(2) - ppiclf_binb(1))/d2new(1))
      ppiclf_ndygp = floor( (ppiclf_binb(4) - ppiclf_binb(3))/d2new(2))
      ppiclf_ndzgp = 1
      if (ppiclf_rparam(12) .gt. 2) ppiclf_ndzgp = 
     >                floor( (ppiclf_binb(6) - ppiclf_binb(5))/d2new(3))


      ! dz comment 3/9/2019
c     if (ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp .gt. ppiclf_np .or. 
c    >    int(ppiclf_rparam(4)) .eq. 1) then
         nmax = 1000
         d2chk_save = ppiclf_d2chk(2)
         
         do i=1,nmax
         do j=0,ndim-1
            ifac(j+1) = 1 + i
            d2new(j+1) = (ppiclf_binb(2+2*j) - ppiclf_binb(1+2*j))/
     >                    ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)

            if( nbb .gt. ppiclf_np ) then
            ! dz comment 3/9/2019
c           if( int(ppiclf_rparam(4)) .eq. 1 .or.
c    >          int(ppiclf_rparam(4)) .eq. 0 .and.d2new(j+1).lt.d2chk_save)
c    >          then
               icount(j+1) = icount(j+1) + 1
               ifac(j+1) = ifac(j+1) - icount(j+1)
               d2new(j+1) = (ppiclf_binb(2+2*j) -ppiclf_binb(1+2*j))/
     >                       ifac(j+1)
c           endif
            endif
         enddo
            if (icount(1) .gt. 0) then
            if (icount(2) .gt. 0) then
            if (icount(3) .gt. 0) then
               exit
            endif
            endif
            endif
         enddo
c     endif

! -------------------------------------------------------
c SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! how many spacings in each direction
      ppiclf_ndxgp = floor( (ppiclf_binb(2) - ppiclf_binb(1))/d2new(1))
      ppiclf_ndygp = floor( (ppiclf_binb(4) - ppiclf_binb(3))/d2new(2))
      ppiclf_ndzgp = 1
      if (ppiclf_rparam(12) .gt. 2) ppiclf_ndzgp = 
     >                floor( (ppiclf_binb(6) - ppiclf_binb(5))/d2new(3))

      
      isize = 4
      call bcast(ppiclf_ndxgp, isize)
      call bcast(ppiclf_ndygp, isize)
      call bcast(ppiclf_ndzgp, isize)
      call bcast(ppiclf_binb , 6*2*isize)

      ! grid spacing for that many spacings
      ppiclf_rdxgp = (ppiclf_binb(2) -ppiclf_binb(1))/real(ppiclf_ndxgp)
      ppiclf_rdygp = (ppiclf_binb(4) -ppiclf_binb(3))/real(ppiclf_ndygp)
      ppiclf_rdzgp = 1.
      if (ppiclf_rparam(12) .gt. 2) 
     >ppiclf_rdzgp = (ppiclf_binb(6) -ppiclf_binb(5))/real(ppiclf_ndzgp)

c     ninc = 2
      rxlbin = ppiclf_binb(1)
      rxrbin = ppiclf_binb(2)
      rylbin = ppiclf_binb(3)
      ryrbin = ppiclf_binb(4)
      rzlbin = ppiclf_binb(5)
      rzrbin = ppiclf_binb(6)
c     if (iperiodicx .ne. 0) then
c        rxlbin = rxlbin - ninc/2*ppiclf_rdxgp
c        rxrbin = rxrbin + ninc/2*ppiclf_rdxgp
c        rxlbin = max(rxlbin,ppiclf_xdrange(1,1))
c        rxrbin = min(rxrbin,ppiclf_xdrange(2,1))
c     endif
c     if (iperiodicy .ne. 0) then
c        rylbin = rylbin - ninc/2*ppiclf_rdygp
c        ryrbin = ryrbin + ninc/2*ppiclf_rdygp
c        rylbin = max(rylbin,ppiclf_xdrange(1,2))
c        ryrbin = min(ryrbin,ppiclf_xdrange(2,2))
c     endif
c     if (iperiodicz .ne. 0) then
c     if (ppiclf_rparam(12) .gt. 2) then
c        rzlbin = rzlbin - ninc/2*ppiclf_rdzgp
c        rzrbin = rzrbin + ninc/2*ppiclf_rdzgp
c        rzlbin = max(rzlbin,ppiclf_xdrange(1,3))
c        rzrbin = min(rzrbin,ppiclf_xdrange(2,3))
c     endif
c     endif

      ppiclf_nbin = ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp

c     current box coordinates
      if (ppiclf_nid .le. ppiclf_nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_ndxgp)
         jdum = modulo(ppiclf_nid/ppiclf_ndxgp,ppiclf_ndygp)
         kdum = ppiclf_nid/(ppiclf_ndxgp*ppiclf_ndygp)
         ppiclf_binx(1,1) = rxlbin + idum    *ppiclf_rdxgp
         ppiclf_binx(2,1) = rxlbin + (idum+1)*ppiclf_rdxgp
         ppiclf_biny(1,1) = rylbin + jdum    *ppiclf_rdygp
         ppiclf_biny(2,1) = rylbin + (jdum+1)*ppiclf_rdygp
         ppiclf_binz(1,1) = rzlbin + kdum    *ppiclf_rdzgp
         ppiclf_binz(2,1) = rzlbin + (kdum+1)*ppiclf_rdzgp
    
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_rdxgp/ppiclf_rparam(5)) + 1 + 1
         ppiclf_by = floor(ppiclf_rdygp/ppiclf_rparam(5)) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_rparam(12) .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_rdzgp/ppiclf_rparam(5)) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_rparam(6)
         ppiclf_by = ppiclf_by*ppiclf_rparam(6)
         if (ppiclf_rparam(12) .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_rparam(6)

c        ! force for now!!! remove and uncomment above
c        ppiclf_bx = 3
c        ppiclf_by = 3
c        ppiclf_bz = 3

         if ((ppiclf_bx .gt. PPICLF_BX1) .or. 
     >       (ppiclf_by .gt. PPICLF_BY1) .or.
     >       (ppiclf_bz .gt. PPICLF_BZ1)) then
               if (ppiclf_nid .eq. 0) write(6, *) 
     >          'INCREASE GRID ALLOCATION PPICLF_B*1'
     >         , ppiclf_bx,ppiclf_by,ppiclf_bz
               return
         endif

         ppiclf_rdx = ppiclf_rdxgp/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_rdygp/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_rparam(12) .gt. 2) 
     >      ppiclf_rdz = ppiclf_rdzgp/(ppiclf_bz-1)

         ndumx = ppiclf_ndxgp*(ppiclf_bx-1) + 1
         ndumy = ppiclf_ndygp*(ppiclf_by-1) + 1
         ndumz = ppiclf_ndzgp*(ppiclf_bz-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = ppiclf_binx(1,1) + (i-1)*ppiclf_rdx
            ppiclf_grid_y(i,j,k) = ppiclf_biny(1,1) + (j-1)*ppiclf_rdy
            ppiclf_grid_z(i,j,k) = ppiclf_binz(1,1) + (k-1)*ppiclf_rdz

            ! indicies global, note new mapping
c           ppiclf_grid_i(i,j,k) = ppiclf_nid*ppiclf_bx*ppiclf_by*ppiclf_bz +
c    >                        (i-1) + ppiclf_bx*(j-1) + ppiclf_bx*ppiclf_by*(k-1)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp
            ppiclf_grid_ii(i,j,k) = itmp
            ppiclf_grid_jj(i,j,k) = jtmp
            ppiclf_grid_kk(i,j,k) = ktmp

         enddo
         enddo
         enddo

      endif

c     rmax = glmax(ppiclf_binz(1,1),2)
c     if (ppiclf_nid .eq. 0) write(6,*) 'MAXX:', rmax

c     if (int(ppiclf_rparam(4)) .eq. 1) return ! only for projection

! see which bins are in which elements
c     ppiclf_neltb = 0
c     do ie=1,ppiclf_nelt
c     do k=1,ppiclf_lz1
c     do j=1,ppiclf_ly1
c     do i=1,ppiclf_lx1
c        rxval = xm1(i,j,k,ie) 
c        ryval = ym1(i,j,k,ie) 
c        rzval = 0.
c        if(if3d) rzval = zm1(i,j,k,ie)

c        if (rxval .gt. ppiclf_binb(2)) goto 1233
c        if (rxval .lt. ppiclf_binb(1)) goto 1233
c        if (ryval .gt. ppiclf_binb(4)) goto 1233
c        if (ryval .lt. ppiclf_binb(3)) goto 1233
c        if (if3d .and. rzval .gt. ppiclf_binb(6)) goto 1233
c        if (if3d .and. rzval .lt. ppiclf_binb(5)) goto 1233

c        ii    = floor((rxval-ppiclf_binb(1))/ppiclf_rdxgp) 
c        jj    = floor((ryval-ppiclf_binb(3))/ppiclf_rdygp) 
c        kk    = floor((rzval-ppiclf_binb(5))/ppiclf_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
c        if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
c        if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
c        ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk
c        nrank = modulo(ndum,np)

c        ppiclf_neltb = ppiclf_neltb + 1
c        if(ppiclf_neltb .gt. ppiclf_lbmax) then
c          write(6,*) 'increase lbmax',nid,ppiclf_neltb,ppiclf_lbmax
c          call exitt
c        endif

c        ppiclf_er_map(1,ppiclf_neltb) = ie
c        ppiclf_er_map(2,ppiclf_neltb) = nid
c        ppiclf_er_map(3,ppiclf_neltb) = ndum
c        ppiclf_er_map(4,ppiclf_neltb) = nrank
c        ppiclf_er_map(5,ppiclf_neltb) = nrank
c        ppiclf_er_map(6,ppiclf_neltb) = nrank

c        if (ppiclf_neltb .gt. 1) then
c        do il=1,ppiclf_neltb-1
c           if (ppiclf_er_map(1,il) .eq. ie) then
c           if (ppiclf_er_map(4,il) .eq. nrank) then
c              ppiclf_neltb = ppiclf_neltb - 1
c              goto 1233
c           endif
c           endif
c        enddo
c        endif
c1233 continue
c     enddo
c     enddo
c     enddo
c     enddo

c     nxyz = lx1*ly1*lz1
c     do ie=1,ppiclf_neltb
c        iee = ppiclf_er_map(1,ie)
c        call copy(ppiclf_xm1b(1,1,1,1,ie), xm1(1,1,1,iee),nxyz)
c        call copy(ppiclf_xm1b(1,1,1,2,ie), ym1(1,1,1,iee),nxyz)
c        call copy(ppiclf_xm1b(1,1,1,3,ie), zm1(1,1,1,iee),nxyz)
c     enddo

c     ppiclf_neltbb = ppiclf_neltb

c     do ie=1,ppiclf_neltbb
c        call icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie),PPICLF_LRMAX)
c     enddo


c     nl   = 0
c     nii  = PPICLF_LRMAX
c     njj  = 6
c     nrr  = nxyz*3
c     nkey = 3
c     call fgslib_crystal_tuple_transfer(i_cr_hndl,ppiclf_neltb,ppiclf_lbmax
c    >                  , ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
c     call fgslib_crystal_tuple_sort    (i_cr_hndl,ppiclf_neltb
c    $              , ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,1)


c     do ie=1,ppiclf_neltb
c     do k=1,nz1
c     do j=1,ny1
c     do i=1,nx1
c        rxval = ppiclf_xm1b(i,j,k,1,ie)
c        ryval = ppiclf_xm1b(i,j,k,2,ie)
c        rzval = 0.
c        if(if3d) rzval = ppiclf_xm1b(i,j,k,3,ie)
c        
c        ii    = floor((rxval-ppiclf_binb(1))/ppiclf_rdxgp) 
c        jj    = floor((ryval-ppiclf_binb(3))/ppiclf_rdygp) 
c        kk    = floor((rzval-ppiclf_binb(5))/ppiclf_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
c        if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
c        if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
c        ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk

c        ppiclf_modgp(i,j,k,ie,1) = ii
c        ppiclf_modgp(i,j,k,ie,2) = jj
c        ppiclf_modgp(i,j,k,ie,3) = kk
c        ppiclf_modgp(i,j,k,ie,4) = ndum
c  
c     enddo
c     enddo
c     enddo
c     enddo

c     do ie=1,ppiclf_neltb
c        ppiclf_xerange(1,1,ie) = vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
c        ppiclf_xerange(2,1,ie) = vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
c        ppiclf_xerange(1,2,ie) = vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
c        ppiclf_xerange(2,2,ie) = vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
c        ppiclf_xerange(1,3,ie) = vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
c        ppiclf_xerange(2,3,ie) = vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)

c        ilow  = floor((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/ppiclf_rdxgp)
c        ihigh = floor((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/ppiclf_rdxgp)
c        jlow  = floor((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/ppiclf_rdygp)
c        jhigh = floor((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/ppiclf_rdygp)
c        klow  = floor((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/ppiclf_rdzgp)
c        khigh = floor((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/ppiclf_rdzgp)
c        if (.not. if3d) then
c           klow = 0
c           khigh = 0
c        endif

c        ppiclf_el_map(1,ie) = ilow  + ppiclf_ndxgp*jlow  
c    >                            + ppiclf_ndxgp*ppiclf_ndygp*klow
c        ppiclf_el_map(2,ie) = ihigh + ppiclf_ndxgp*jhigh 
c    >                            + ppiclf_ndxgp*ppiclf_ndygp*khigh
c        ppiclf_el_map(3,ie) = ilow
c        ppiclf_el_map(4,ie) = ihigh
c        ppiclf_el_map(5,ie) = jlow
c        ppiclf_el_map(6,ie) = jhigh
c        ppiclf_el_map(7,ie) = klow
c        ppiclf_el_map(8,ie) = khigh
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_findpts
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

c     common /intp_h/ ih_intp(2,1)

c     ih_intp2 = ih_intp(2,i_fp_hndl)

      ix = 1
      iy = 2
      iz = 3

c     call fgslib_findpts(ih_intp2           !   call fgslib_findpts( ihndl,
c    $        , ppiclf_iprop (1 ,1),PPICLF_LIP        !   $             rcode,1,
c    $        , ppiclf_iprop (3 ,1),PPICLF_LIP        !   &             proc,1,
c    $        , ppiclf_iprop (2 ,1),PPICLF_LIP        !   &             elid,1,
c    $        , ppiclf_rprop2(1 ,1),PPICLF_LRP2       !   &             rst,ndim,
c    $        , ppiclf_rprop2(4 ,1),PPICLF_LRP2       !   &             dist,1,
c    $        , ppiclf_y     (ix,1),PPICLF_LRS        !   &             pts(    1),1,
c    $        , ppiclf_y     (iy,1),PPICLF_LRS        !   &             pts(  n+1),1,
c    $        , ppiclf_y     (iz,1),PPICLF_LRS ,PPICLF_NPART) !   &             pts(2*n+1),1,n)

      do i=1,ppiclf_npart
         ppiclf_iprop(4,i) = ppiclf_iprop(3,i)
      enddo

      ! instead get which bin it is in
      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ! add below....

         ii    = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_rdxgp) 
         jj    = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_rdygp) 
         kk    = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_rdzgp) 
         if (ppiclf_rparam(12) .lt. 3) kk = 0
         if (ii .eq. ppiclf_ndxgp) ii = ppiclf_ndxgp - 1
         if (jj .eq. ppiclf_ndygp) jj = ppiclf_ndygp - 1
         if (kk .eq. ppiclf_ndzgp) kk = ppiclf_ndzgp - 1
         ndum  = ii + ppiclf_ndxgp*jj + ppiclf_ndxgp*ppiclf_ndygp*kk
         nrank = modulo(ndum, ppiclf_np)

         ppiclf_iprop(8,i)  = ii
         ppiclf_iprop(9,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
      enddo


      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_crystal
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      logical partl    
      integer ppiclf_ipmap1(1,PPICLF_LPART)
     >       ,ppiclf_ipmap2(1,PPICLF_LPART)
     >       ,ppiclf_ipmap3(1,PPICLF_LPART)
     >       ,ppiclf_ipmap4(1,PPICLF_LPART)

      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2)
      real rwork(lrf,PPICLF_LPART)

      do i=1,ppiclf_npart
         ic = 1
         call copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
      enddo

      j0 = 4
      call fgslib_crystal_tuple_transfer(i_cr_hndl,ppiclf_npart 
     >      ,PPICLF_LPART,ppiclf_iprop ,PPICLF_LIP,partl,0,rwork,lrf,j0)

      do i=1,ppiclf_npart
         ic = 1
         call copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
      enddo
        
      return
      end
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_ghost_create
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      character*132 deathmessage
      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),gpsave(27)
      real map(PPICLF_LRP_PRO)

      integer  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp

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

      if (ppiclf_rparam(12) .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = int(ppiclf_rparam(8))
      iperiodicy = int(ppiclf_rparam(9))
      iperiodicz = int(ppiclf_rparam(10))

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
      zdlen = -1.
      if (ppiclf_rparam(12) .gt. 2) 
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

         call ppiclf_project_map(map,ppiclf_y(1,ip),ppiclf_ydot(1,ip)
     >                       ,ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))
         ppiclf_cp_map(1,ip) = ppiclf_y(jx,ip)      ! x coord
         ppiclf_cp_map(2,ip) = ppiclf_y(jy,ip)      ! y coord
         ppiclf_cp_map(3,ip) = ppiclf_y(jz,ip)      ! z coord
         do j=1,PPICLF_LRP_PRO
            ppiclf_cp_map(3+j,ip) = map(j)
         enddo

         rxval = ppiclf_cp_map(1,ip)
         ryval = ppiclf_cp_map(2,ip)
         rzval = 0.
         if (ppiclf_rparam(12) .gt. 2) rzval = ppiclf_cp_map(3,ip)

         iip    = ppiclf_iprop(8,ip)
         jjp    = ppiclf_iprop(9,ip)
         kkp    = ppiclf_iprop(10,ip)

         rxl = ppiclf_binb(1) + ppiclf_rdxgp*iip
         rxr = rxl + ppiclf_rdxgp
         ryl = ppiclf_binb(3) + ppiclf_rdygp*jjp
         ryr = ryl + ppiclf_rdygp
         rzl = 0.0
         rzr = 0.0
         if (ppiclf_rparam(12) .gt. 2) then
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
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
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

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

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

            call ppiclf_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
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
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
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

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

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

            call ppiclf_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
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
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(2))**2
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

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

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

            call ppiclf_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
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
      subroutine ppiclf_comm_check_periodic_gp(rxnew,rxdrng,iadd)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
c
      real rxnew(3), rxdrng(3)
      integer iadd(3), irett(3), ntype, ntypel(7)

      xloc = rxnew(1)
      yloc = rxnew(2)
      zloc = rxnew(3)

      xdlen = rxdrng(1)
      ydlen = rxdrng(2)
      zdlen = rxdrng(3)

      ii = iadd(1)
      jj = iadd(2)
      kk = iadd(3)

      irett(1) = 0
      irett(2) = 0
      irett(3) = 0

      if (xdlen .gt. 0 ) then
      if (ii .ge. ppiclf_ndxgp) then
         xloc = xloc - xdlen
         irett(1) = 1
         goto 123
      endif
      endif
      if (xdlen .gt. 0 ) then
      if (ii .lt. 0) then
         xloc = xloc + xdlen
         irett(1) = 1
         goto 123
      endif
      endif

  123 continue    
      if (ydlen .gt. 0 ) then
      if (jj .ge. ppiclf_ndygp) then
         yloc = yloc - ydlen
         irett(2) = 1
         goto 124
      endif
      endif
      if (ydlen .gt. 0 ) then
      if (jj .lt. 0) then
         yloc = yloc + ydlen
         irett(2) = 1
         goto 124
      endif
      endif
  124 continue

      if (ppiclf_rparam(12) .gt. 2) then
         if (zdlen .gt. 0 ) then
         if (kk .ge. ppiclf_ndzgp) then
            zloc = zloc - zdlen
            irett(3) = 1
            goto 125
         endif
         endif
         if (zdlen .gt. 0 ) then
         if (kk .lt. 0) then
            zloc = zloc + zdlen
            irett(3) = 1
            goto 125
         endif
         endif
      endif
  125 continue

      rxnew(1) = xloc
      rxnew(2) = yloc
      rxnew(3) = zloc

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_ghost_send
c
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"

      logical partl         

      ! send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
