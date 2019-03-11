!-----------------------------------------------------------------------
      subroutine lpm_comm_setup
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

c     call lpm_comm_interp_setup(i_fp_hndl,0.0,idum,lpm_nelt)
      call fgslib_crystal_setup(i_cr_hndl,lpm_comm,lpm_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_comm_interp_setup(ih,tolin,nmsh,nelm)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      ! need to add xm1, etc in LPM for this case

c     common /intp_h/ ih_intp(2,1)
c     common /intp/   tol

c     data ihcounter /0/
c     save ihcounter

c     real xmi, ymi, zmi
c     common /CBXMI/ xmi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT),
c    $               ymi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT),
c    $               zmi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT)
c
c     real w(2*LPM_LX1*LPM_LY1*LPM_LZ1)
c     logical if3d

c     tol = max(5e-13,tolin)
c     npt_max = 128
c     bb_t    = 0.01

c     ! fix below
c     if3d = .true.

c     if (nmsh.ge.1 .and. nmsh.lt.LPM_LX1-1) then
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
      subroutine lpm_comm_bin_setup
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      integer  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer  ifac(3), icount(3)
      real     d2new(3)
      logical  partl

      real                      lpm_xerange(2,3,lpm_lbmax)
      common /lpm_elementrange/ lpm_xerange

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

      if (lpm_rparam(12) .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 3

      iperiodicx = int(lpm_rparam(8))
      iperiodicy = int(lpm_rparam(9))
      iperiodicz = int(lpm_rparam(10))
      ndim       = int(lpm_rparam(12))

      ! compute binb
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      do i=1,lpm_npart
         rduml = lpm_y(ix,i) - lpm_d2chk(2)
         rdumr = lpm_y(ix,i) + lpm_d2chk(2)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = lpm_y(iy,i) - lpm_d2chk(2)
         rdumr = lpm_y(iy,i) + lpm_d2chk(2)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         rduml = lpm_y(iz,i) - lpm_d2chk(2)
         rdumr = lpm_y(iz,i) + lpm_d2chk(2)
         if (rduml .lt. zmin) zmin = rduml
         if (rdumr .gt. zmax) zmax = rdumr
      enddo

      lpm_binb(1) = glmin(xmin,1)
      lpm_binb(2) = glmax(xmax,1)
      lpm_binb(3) = glmin(ymin,1)
      lpm_binb(4) = glmax(ymax,1)
      lpm_binb(5) = 0.0
      lpm_binb(6) = 0.0
      if(lpm_rparam(12) .gt. 2) lpm_binb(5) = glmin(zmin,1)
      if(lpm_rparam(12) .gt. 2) lpm_binb(6) = glmax(zmax,1)


      ! comment for now.....
c     lpm_binb(1) = max(lpm_binb(1),lpm_xdrange(1,1))
c     lpm_binb(2) = min(lpm_binb(2),lpm_xdrange(2,1))
c     lpm_binb(3) = max(lpm_binb(3),lpm_xdrange(1,2))
c     lpm_binb(4) = min(lpm_binb(4),lpm_xdrange(2,2))
c     if(lpm_rparam(12) .gt. 2) lpm_binb(5) = 
c    >                           max(lpm_binb(5),lpm_xdrange(1,3))
c     if(lpm_rparam(12) .gt. 2) lpm_binb(6) = 
c    >                           min(lpm_binb(6),lpm_xdrange(2,3))
c     if (iperiodicx .eq. 0) then
c        lpm_binb(1) = lpm_xdrange(1,1)
c        lpm_binb(2) = lpm_xdrange(2,1)
c     endif
c     if (iperiodicy .eq. 0) then
c        lpm_binb(3) = lpm_xdrange(1,2)
c        lpm_binb(4) = lpm_xdrange(2,2)
c     endif
c     if (iperiodicz .eq. 0 .and. lpm_rparam(12) .gt. 2) then
c        lpm_binb(5) = lpm_xdrange(1,3)
c        lpm_binb(6) = lpm_xdrange(2,3)
c     endif

      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      icount(1) = 0
      icount(2) = 0
      icount(3) = 0
      d2new(1) = lpm_d2chk(2)
      d2new(2) = lpm_d2chk(2)
      d2new(3) = lpm_d2chk(2)

      lpm_ndxgp = floor( (lpm_binb(2) - lpm_binb(1))/d2new(1))
      lpm_ndygp = floor( (lpm_binb(4) - lpm_binb(3))/d2new(2))
      lpm_ndzgp = 1
      if (lpm_rparam(12) .gt. 2) lpm_ndzgp = 
     >                      floor( (lpm_binb(6) - lpm_binb(5))/d2new(3))


      ! dz comment 3/9/2019
c     if (lpm_ndxgp*lpm_ndygp*lpm_ndzgp .gt. lpm_np .or. 
c    >    int(lpm_rparam(4)) .eq. 1) then
         nmax = 1000
         d2chk_save = lpm_d2chk(2)
         
         do i=1,nmax
         do j=0,ndim-1
            ifac(j+1) = 1 + i
            d2new(j+1) = (lpm_binb(2+2*j) - lpm_binb(1+2*j))/ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)

            if( nbb .gt. lpm_np ) then
            ! dz comment 3/9/2019
c           if( int(lpm_rparam(4)) .eq. 1 .or.
c    >          int(lpm_rparam(4)) .eq. 0 .and.d2new(j+1).lt.d2chk_save)
c    >          then
               icount(j+1) = icount(j+1) + 1
               ifac(j+1) = ifac(j+1) - icount(j+1)
               d2new(j+1) = (lpm_binb(2+2*j) -lpm_binb(1+2*j))/ifac(j+1)
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
      lpm_ndxgp = floor( (lpm_binb(2) - lpm_binb(1))/d2new(1))
      lpm_ndygp = floor( (lpm_binb(4) - lpm_binb(3))/d2new(2))
      lpm_ndzgp = 1
      if (lpm_rparam(12) .gt. 2) lpm_ndzgp = 
     >                      floor( (lpm_binb(6) - lpm_binb(5))/d2new(3))

      
      isize = 4
      call bcast(lpm_ndxgp, isize)
      call bcast(lpm_ndygp, isize)
      call bcast(lpm_ndzgp, isize)
      call bcast(lpm_binb , 6*2*isize)

      ! grid spacing for that many spacings
      lpm_rdxgp = (lpm_binb(2) - lpm_binb(1))/real(lpm_ndxgp)
      lpm_rdygp = (lpm_binb(4) - lpm_binb(3))/real(lpm_ndygp)
      lpm_rdzgp = 1.
      if (lpm_rparam(12) .gt. 2) lpm_rdzgp = (lpm_binb(6) - lpm_binb(5))
     >                           /real(lpm_ndzgp)

c     ninc = 2
      rxlbin = lpm_binb(1)
      rxrbin = lpm_binb(2)
      rylbin = lpm_binb(3)
      ryrbin = lpm_binb(4)
      rzlbin = lpm_binb(5)
      rzrbin = lpm_binb(6)
c     if (iperiodicx .ne. 0) then
c        rxlbin = rxlbin - ninc/2*lpm_rdxgp
c        rxrbin = rxrbin + ninc/2*lpm_rdxgp
c        rxlbin = max(rxlbin,lpm_xdrange(1,1))
c        rxrbin = min(rxrbin,lpm_xdrange(2,1))
c     endif
c     if (iperiodicy .ne. 0) then
c        rylbin = rylbin - ninc/2*lpm_rdygp
c        ryrbin = ryrbin + ninc/2*lpm_rdygp
c        rylbin = max(rylbin,lpm_xdrange(1,2))
c        ryrbin = min(ryrbin,lpm_xdrange(2,2))
c     endif
c     if (iperiodicz .ne. 0) then
c     if (lpm_rparam(12) .gt. 2) then
c        rzlbin = rzlbin - ninc/2*lpm_rdzgp
c        rzrbin = rzrbin + ninc/2*lpm_rdzgp
c        rzlbin = max(rzlbin,lpm_xdrange(1,3))
c        rzrbin = min(rzrbin,lpm_xdrange(2,3))
c     endif
c     endif

      lpm_nbin = lpm_ndxgp*lpm_ndygp*lpm_ndzgp

c     current box coordinates
      if (lpm_nid .le. lpm_nbin-1) then
         idum = modulo(lpm_nid,lpm_ndxgp)
         jdum = modulo(lpm_nid/lpm_ndxgp,lpm_ndygp)
         kdum = lpm_nid/(lpm_ndxgp*lpm_ndygp)
         lpm_binx(1,1) = rxlbin + idum    *lpm_rdxgp
         lpm_binx(2,1) = rxlbin + (idum+1)*lpm_rdxgp
         lpm_biny(1,1) = rylbin + jdum    *lpm_rdygp
         lpm_biny(2,1) = rylbin + (jdum+1)*lpm_rdygp
         lpm_binz(1,1) = rzlbin + kdum    *lpm_rdzgp
         lpm_binz(2,1) = rzlbin + (kdum+1)*lpm_rdzgp
    
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         lpm_bx = floor(lpm_rdxgp/lpm_rparam(5)) + 1 + 1
         lpm_by = floor(lpm_rdygp/lpm_rparam(5)) + 1 + 1
         lpm_bz = 1
         if (lpm_rparam(12) .gt. 2) 
     >      lpm_bz = floor(lpm_rdzgp/lpm_rparam(5)) + 1 + 1

         lpm_bx = lpm_bx*lpm_rparam(6)
         lpm_by = lpm_by*lpm_rparam(6)
         if (lpm_rparam(12) .gt. 2) 
     >      lpm_bz = lpm_bz*lpm_rparam(6)

c        ! force for now!!! remove and uncomment above
c        lpm_bx = 3
c        lpm_by = 3
c        lpm_bz = 3

         if ((lpm_bx .gt. LPM_BX1) .or. 
     >       (lpm_by .gt. LPM_BY1) .or.
     >       (lpm_bz .gt. LPM_BZ1)) then
               if (lpm_nid .eq. 0) write(6, *) 
     >          'INCREASE GRID ALLOCATION LPM_B*1', lpm_bx,lpm_by,lpm_bz
               return
         endif

         lpm_rdx = lpm_rdxgp/(lpm_bx-1)
         lpm_rdy = lpm_rdygp/(lpm_by-1)
         lpm_rdz = 0
         if (lpm_rparam(12) .gt. 2) 
     >      lpm_rdz = lpm_rdzgp/(lpm_bz-1)

         ndumx = lpm_ndxgp*(lpm_bx-1) + 1
         ndumy = lpm_ndygp*(lpm_by-1) + 1
         ndumz = lpm_ndzgp*(lpm_bz-1) + 1
    
         do k=1,lpm_bz
         do j=1,lpm_by
         do i=1,lpm_bx
            lpm_grid_x(i,j,k) = lpm_binx(1,1) + (i-1)*lpm_rdx
            lpm_grid_y(i,j,k) = lpm_biny(1,1) + (j-1)*lpm_rdy
            lpm_grid_z(i,j,k) = lpm_binz(1,1) + (k-1)*lpm_rdz

            ! indicies global, note new mapping
c           lpm_grid_i(i,j,k) = lpm_nid*lpm_bx*lpm_by*lpm_bz +
c    >                        (i-1) + lpm_bx*(j-1) + lpm_bx*lpm_by*(k-1)

            itmp = idum*(lpm_bx-1) + (i-1)
            jtmp = jdum*(lpm_by-1) + (j-1)
            ktmp = kdum*(lpm_bz-1) + (k-1)
    
            lpm_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp
            lpm_grid_ii(i,j,k) = itmp
            lpm_grid_jj(i,j,k) = jtmp
            lpm_grid_kk(i,j,k) = ktmp

         enddo
         enddo
         enddo

      endif

c     rmax = glmax(lpm_binz(1,1),2)
c     if (lpm_nid .eq. 0) write(6,*) 'MAXX:', rmax

c     if (int(lpm_rparam(4)) .eq. 1) return ! only for projection

! see which bins are in which elements
c     lpm_neltb = 0
c     do ie=1,lpm_nelt
c     do k=1,lpm_lz1
c     do j=1,lpm_ly1
c     do i=1,lpm_lx1
c        rxval = xm1(i,j,k,ie) 
c        ryval = ym1(i,j,k,ie) 
c        rzval = 0.
c        if(if3d) rzval = zm1(i,j,k,ie)

c        if (rxval .gt. lpm_binb(2)) goto 1233
c        if (rxval .lt. lpm_binb(1)) goto 1233
c        if (ryval .gt. lpm_binb(4)) goto 1233
c        if (ryval .lt. lpm_binb(3)) goto 1233
c        if (if3d .and. rzval .gt. lpm_binb(6)) goto 1233
c        if (if3d .and. rzval .lt. lpm_binb(5)) goto 1233

c        ii    = floor((rxval-lpm_binb(1))/lpm_rdxgp) 
c        jj    = floor((ryval-lpm_binb(3))/lpm_rdygp) 
c        kk    = floor((rzval-lpm_binb(5))/lpm_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
c        if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
c        if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
c        ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk
c        nrank = modulo(ndum,np)

c        lpm_neltb = lpm_neltb + 1
c        if(lpm_neltb .gt. lpm_lbmax) then
c          write(6,*) 'increase lbmax',nid,lpm_neltb,lpm_lbmax
c          call exitt
c        endif

c        lpm_er_map(1,lpm_neltb) = ie
c        lpm_er_map(2,lpm_neltb) = nid
c        lpm_er_map(3,lpm_neltb) = ndum
c        lpm_er_map(4,lpm_neltb) = nrank
c        lpm_er_map(5,lpm_neltb) = nrank
c        lpm_er_map(6,lpm_neltb) = nrank

c        if (lpm_neltb .gt. 1) then
c        do il=1,lpm_neltb-1
c           if (lpm_er_map(1,il) .eq. ie) then
c           if (lpm_er_map(4,il) .eq. nrank) then
c              lpm_neltb = lpm_neltb - 1
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
c     do ie=1,lpm_neltb
c        iee = lpm_er_map(1,ie)
c        call copy(lpm_xm1b(1,1,1,1,ie), xm1(1,1,1,iee),nxyz)
c        call copy(lpm_xm1b(1,1,1,2,ie), ym1(1,1,1,iee),nxyz)
c        call copy(lpm_xm1b(1,1,1,3,ie), zm1(1,1,1,iee),nxyz)
c     enddo

c     lpm_neltbb = lpm_neltb

c     do ie=1,lpm_neltbb
c        call icopy(lpm_er_maps(1,ie),lpm_er_map(1,ie),LPM_LRMAX)
c     enddo


c     nl   = 0
c     nii  = LPM_LRMAX
c     njj  = 6
c     nrr  = nxyz*3
c     nkey = 3
c     call fgslib_crystal_tuple_transfer(i_cr_hndl,lpm_neltb,lpm_lbmax
c    >                  , lpm_er_map,nii,partl,nl,lpm_xm1b,nrr,njj)
c     call fgslib_crystal_tuple_sort    (i_cr_hndl,lpm_neltb
c    $              , lpm_er_map,nii,partl,nl,lpm_xm1b,nrr,nkey,1)


c     do ie=1,lpm_neltb
c     do k=1,nz1
c     do j=1,ny1
c     do i=1,nx1
c        rxval = lpm_xm1b(i,j,k,1,ie)
c        ryval = lpm_xm1b(i,j,k,2,ie)
c        rzval = 0.
c        if(if3d) rzval = lpm_xm1b(i,j,k,3,ie)
c        
c        ii    = floor((rxval-lpm_binb(1))/lpm_rdxgp) 
c        jj    = floor((ryval-lpm_binb(3))/lpm_rdygp) 
c        kk    = floor((rzval-lpm_binb(5))/lpm_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
c        if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
c        if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
c        ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk

c        lpm_modgp(i,j,k,ie,1) = ii
c        lpm_modgp(i,j,k,ie,2) = jj
c        lpm_modgp(i,j,k,ie,3) = kk
c        lpm_modgp(i,j,k,ie,4) = ndum
c  
c     enddo
c     enddo
c     enddo
c     enddo

c     do ie=1,lpm_neltb
c        lpm_xerange(1,1,ie) = vlmin(lpm_xm1b(1,1,1,1,ie),nxyz)
c        lpm_xerange(2,1,ie) = vlmax(lpm_xm1b(1,1,1,1,ie),nxyz)
c        lpm_xerange(1,2,ie) = vlmin(lpm_xm1b(1,1,1,2,ie),nxyz)
c        lpm_xerange(2,2,ie) = vlmax(lpm_xm1b(1,1,1,2,ie),nxyz)
c        lpm_xerange(1,3,ie) = vlmin(lpm_xm1b(1,1,1,3,ie),nxyz)
c        lpm_xerange(2,3,ie) = vlmax(lpm_xm1b(1,1,1,3,ie),nxyz)

c        ilow  = floor((lpm_xerange(1,1,ie) - lpm_binb(1))/lpm_rdxgp)
c        ihigh = floor((lpm_xerange(2,1,ie) - lpm_binb(1))/lpm_rdxgp)
c        jlow  = floor((lpm_xerange(1,2,ie) - lpm_binb(3))/lpm_rdygp)
c        jhigh = floor((lpm_xerange(2,2,ie) - lpm_binb(3))/lpm_rdygp)
c        klow  = floor((lpm_xerange(1,3,ie) - lpm_binb(5))/lpm_rdzgp)
c        khigh = floor((lpm_xerange(2,3,ie) - lpm_binb(5))/lpm_rdzgp)
c        if (.not. if3d) then
c           klow = 0
c           khigh = 0
c        endif

c        lpm_el_map(1,ie) = ilow  + lpm_ndxgp*jlow  
c    >                            + lpm_ndxgp*lpm_ndygp*klow
c        lpm_el_map(2,ie) = ihigh + lpm_ndxgp*jhigh 
c    >                            + lpm_ndxgp*lpm_ndygp*khigh
c        lpm_el_map(3,ie) = ilow
c        lpm_el_map(4,ie) = ihigh
c        lpm_el_map(5,ie) = jlow
c        lpm_el_map(6,ie) = jhigh
c        lpm_el_map(7,ie) = klow
c        lpm_el_map(8,ie) = khigh
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_comm_findpts
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

c     common /intp_h/ ih_intp(2,1)

c     ih_intp2 = ih_intp(2,i_fp_hndl)

      ix = 1
      iy = 2
      iz = 3

c     call fgslib_findpts(ih_intp2           !   call fgslib_findpts( ihndl,
c    $        , lpm_iprop (1 ,1),LPM_LIP        !   $             rcode,1,
c    $        , lpm_iprop (3 ,1),LPM_LIP        !   &             proc,1,
c    $        , lpm_iprop (2 ,1),LPM_LIP        !   &             elid,1,
c    $        , lpm_rprop2(1 ,1),LPM_LRP2       !   &             rst,ndim,
c    $        , lpm_rprop2(4 ,1),LPM_LRP2       !   &             dist,1,
c    $        , lpm_y     (ix,1),LPM_LRS        !   &             pts(    1),1,
c    $        , lpm_y     (iy,1),LPM_LRS        !   &             pts(  n+1),1,
c    $        , lpm_y     (iz,1),LPM_LRS ,LPM_NPART) !   &             pts(2*n+1),1,n)

      do i=1,lpm_npart
         lpm_iprop(4,i) = lpm_iprop(3,i)
      enddo

      ! instead get which bin it is in
      do i=1,lpm_npart
         ! check if particles are greater or less than binb bounds....
         ! add below....

         ii    = floor((lpm_y(ix,i)-lpm_binb(1))/lpm_rdxgp) 
         jj    = floor((lpm_y(iy,i)-lpm_binb(3))/lpm_rdygp) 
         kk    = floor((lpm_y(iz,i)-lpm_binb(5))/lpm_rdzgp) 
         if (lpm_rparam(12) .lt. 3) kk = 0
         if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
         if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
         if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
         ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk
         nrank = modulo(ndum, lpm_np)

         lpm_iprop(8,i)  = ii
         lpm_iprop(9,i)  = jj
         lpm_iprop(10,i) = kk
         lpm_iprop(11,i) = ndum

         lpm_iprop(4,i)  = nrank ! where particle is actually moved
      enddo


      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_comm_crystal
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      logical partl    
      integer lpm_ipmap1(1,LPM_LPART)
     >       ,lpm_ipmap2(1,LPM_LPART)
     >       ,lpm_ipmap3(1,LPM_LPART)
     >       ,lpm_ipmap4(1,LPM_LPART)

      parameter(lrf = LPM_LRS*4 + LPM_LRP + LPM_LRP2)
      real rwork(lrf,LPM_LPART)

      do i=1,lpm_npart
         ic = 1
         call copy(rwork(ic,i),lpm_y(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_y1((i-1)*LPM_LRS+1),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_ydot(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_ydotc(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_rprop(1,i),LPM_LRP)
         ic = ic + LPM_LRP
         call copy(rwork(ic,i),lpm_rprop2(1,i),LPM_LRP2)
      enddo

      j0 = 4
      call fgslib_crystal_tuple_transfer(i_cr_hndl,lpm_npart ,LPM_LPART
     $           ,lpm_iprop ,LPM_LIP,partl,0,rwork,lrf ,j0)

      do i=1,lpm_npart
         ic = 1
         call copy(lpm_y(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_y1((i-1)*LPM_LRS+1),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_ydot(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_ydotc(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_rprop(1,i),rwork(ic,i),LPM_LRP)
         ic = ic + LPM_LRP
         call copy(lpm_rprop2(1,i),rwork(ic,i),LPM_LRP2)
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
      subroutine lpm_comm_ghost_create
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      character*132 deathmessage
      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),gpsave(27)
      real map(LPM_LRP_PRO)

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

      if (lpm_rparam(12) .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = int(lpm_rparam(8))
      iperiodicy = int(lpm_rparam(9))
      iperiodicz = int(lpm_rparam(10))

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      xdlen = lpm_binb(2) - lpm_binb(1)
      ydlen = lpm_binb(4) - lpm_binb(3)
      zdlen = -1.
      if (lpm_rparam(12) .gt. 2) zdlen = lpm_binb(6) - lpm_binb(5)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      lpm_npart_gp = 0

      rfac = 1.0


      do ip=1,lpm_npart

         call lpm_project_map(map,lpm_y(1,ip),lpm_ydot(1,ip)
     >                       ,lpm_ydotc(1,ip),lpm_rprop(1,ip))
         lpm_cp_map(1,ip) = lpm_y(jx,ip)      ! x coord
         lpm_cp_map(2,ip) = lpm_y(jy,ip)      ! y coord
         lpm_cp_map(3,ip) = lpm_y(jz,ip)      ! z coord
         do j=1,LPM_LRP_PRO
            lpm_cp_map(3+j,ip) = map(j)
         enddo

         rxval = lpm_cp_map(1,ip)
         ryval = lpm_cp_map(2,ip)
         rzval = 0.
         if (lpm_rparam(12) .gt. 2) rzval = lpm_cp_map(3,ip)

         iip    = lpm_iprop(8,ip)
         jjp    = lpm_iprop(9,ip)
         kkp    = lpm_iprop(10,ip)

         rxl = lpm_binb(1) + lpm_rdxgp*iip
         rxr = rxl + lpm_rdxgp
         ryl = lpm_binb(3) + lpm_rdygp*jjp
         ryr = ryl + lpm_rdygp
         rzl = 0.0
         rzr = 0.0
         if (lpm_rparam(12) .gt. 2) then
            rzl = lpm_binb(5) + lpm_rdzgp*kkp
            rzr = rzl + lpm_rdzgp
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
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (lpm_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
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
            if (iig .lt. 0 .or. iig .gt. lpm_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,lpm_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. lpm_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,lpm_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. lpm_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,lpm_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + lpm_ndxgp*jjg + lpm_ndxgp*lpm_ndygp*kkg
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

            call lpm_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
            lpm_npart_gp = lpm_npart_gp + 1
            lpm_iprop_gp(1,lpm_npart_gp) = nrank
            lpm_iprop_gp(2,lpm_npart_gp) = iig
            lpm_iprop_gp(3,lpm_npart_gp) = jjg
            lpm_iprop_gp(4,lpm_npart_gp) = kkg
            lpm_iprop_gp(5,lpm_npart_gp) = ndumn

            lpm_rprop_gp(1,lpm_npart_gp) = rxnew(1)
            lpm_rprop_gp(2,lpm_npart_gp) = rxnew(2)
            lpm_rprop_gp(3,lpm_npart_gp) = rxnew(3)
            do k=4,LPM_LRP_GP
               lpm_rprop_gp(k,lpm_npart_gp) = lpm_cp_map(k,ip)
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
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (lpm_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
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
            if (iig .lt. 0 .or. iig .gt. lpm_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,lpm_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. lpm_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,lpm_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. lpm_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,lpm_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + lpm_ndxgp*jjg + lpm_ndxgp*lpm_ndygp*kkg
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

            call lpm_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
            lpm_npart_gp = lpm_npart_gp + 1
            lpm_iprop_gp(1,lpm_npart_gp) = nrank
            lpm_iprop_gp(2,lpm_npart_gp) = iig
            lpm_iprop_gp(3,lpm_npart_gp) = jjg
            lpm_iprop_gp(4,lpm_npart_gp) = kkg
            lpm_iprop_gp(5,lpm_npart_gp) = ndumn

            lpm_rprop_gp(1,lpm_npart_gp) = rxnew(1)
            lpm_rprop_gp(2,lpm_npart_gp) = rxnew(2)
            lpm_rprop_gp(3,lpm_npart_gp) = rxnew(3)
            do k=4,LPM_LRP_GP
               lpm_rprop_gp(k,lpm_npart_gp) = lpm_cp_map(k,ip)
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
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (lpm_rparam(12) .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*lpm_d2chk(2))**2
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
            if (iig .lt. 0 .or. iig .gt. lpm_ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,lpm_ndxgp)
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. lpm_ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,lpm_ndygp)
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. lpm_ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,lpm_ndzgp)
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + lpm_ndxgp*jjg + lpm_ndxgp*lpm_ndygp*kkg
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

            call lpm_comm_check_periodic_gp(rxnew,rxdrng,iadd)
                 
            lpm_npart_gp = lpm_npart_gp + 1
            lpm_iprop_gp(1,lpm_npart_gp) = nrank
            lpm_iprop_gp(2,lpm_npart_gp) = iig
            lpm_iprop_gp(3,lpm_npart_gp) = jjg
            lpm_iprop_gp(4,lpm_npart_gp) = kkg
            lpm_iprop_gp(5,lpm_npart_gp) = ndumn

            lpm_rprop_gp(1,lpm_npart_gp) = rxnew(1)
            lpm_rprop_gp(2,lpm_npart_gp) = rxnew(2)
            lpm_rprop_gp(3,lpm_npart_gp) = rxnew(3)
            do k=4,LPM_LRP_GP
               lpm_rprop_gp(k,lpm_npart_gp) = lpm_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_comm_check_periodic_gp(rxnew,rxdrng,iadd)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
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
      if (ii .ge. lpm_ndxgp) then
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
      if (jj .ge. lpm_ndygp) then
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

      if (lpm_rparam(12) .gt. 2) then
         if (zdlen .gt. 0 ) then
         if (kk .ge. lpm_ndzgp) then
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
      subroutine lpm_comm_ghost_send
c
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"

      logical partl         

      ! send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl
     >                                  ,lpm_npart_gp,LPM_LPART_GP
     >                                  ,lpm_iprop_gp,LPM_LIP_GP
     >                                  ,partl,0
     >                                  ,lpm_rprop_gp,LPM_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
